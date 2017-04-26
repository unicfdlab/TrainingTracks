/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "weaklyCoupled3DofFsi.H"
#include "volFields.H"
#include "Time.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(weaklyCoupled3DofFsi, 0);
    
    addToRunTimeSelectionTable(functionObject, weaklyCoupled3DofFsi, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::weaklyCoupled3DofFsi::weaklyCoupled3DofFsi
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    forces
    (
        name,
        runTime,
        dict
    ),
    M_(vector::one),
    C_(vector::zero),
    K_(vector::zero),
    R_(1),
    Ymax_(vector::zero),
    append_(false),
    Y_ (vector::zero, vector::zero),
    Yold_(vector::zero, vector::zero),
    coordSys_(NULL)
{
    this->read(dict);
    this->createFsiOutFile(dict);
}

Foam::functionObjects::weaklyCoupled3DofFsi::weaklyCoupled3DofFsi
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    forces
    (
        name,
        obr,
        dict
    ),
    M_(vector::one),
    C_(vector::zero),
    K_(vector::zero),
    R_(1),
    Ymax_(vector::zero),
    append_(false),
    Y_ (vector::zero, vector::zero),
    Yold_(vector::zero, vector::zero),
    coordSys_(NULL)
{
    this->read(dict);
    
    if (Pstream::master())
    {
        List<word> oldFileLines(0);
        if (append_)
        {
            IFstream outOld
            (
                dict.lookup("results")
            );
            
            while(!outOld.eof() && outOld.opened())
            {
                word str(word::null);
                outOld.getLine(str);
                if (!str.empty())
                {
                    oldFileLines.append(str);
                }
            }
        }
        
        this->createFsiOutFile(dict);
        
        if (append_ && oldFileLines.size())
        {
            for(label i=1; i<oldFileLines.size(); i++)
            {
                file(2) << oldFileLines[i] << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::weaklyCoupled3DofFsi::~weaklyCoupled3DofFsi()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::weaklyCoupled3DofFsi::createFsiOutFile(const dictionary& dict)
{
    if (Pstream::master())
    {
        files().resize(3);
        files().set
        (
            2,
            new OFstream
            (
                dict.lookup("results")
            )
        );
        file(2) << "Time;Y1;Y2;Y3;Vy1;Vy2;Vy3;Fy1;Fy2;Fy3" << endl;
    }
}

bool Foam::functionObjects::weaklyCoupled3DofFsi::read(const dictionary& dict)
{
    if (!forces::read(dict))
    {
        return false;
    }
    
    dict.lookup("M") >> M_;
    
    dict.lookup("C") >> C_;
    
    dict.lookup("K") >> K_;
    
    dict.lookup("R") >> R_;
    
    dict.lookup("Ymax") >> Ymax_;
    
    if (coordSys_.empty())
    {
        coordSys_ = coordinateSystem::New (obr_, dict);
    }

    Info << "Reading old state" << endl;
    autoPtr<IOdictionary> weaklyCoupled3DofFsiFsiDictPtr;
    //try to read weaklyCoupled3DofFsiFsi object properties
    {
        volVectorField& yDispl = 
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>("cellDisplacement")
            );
            
            //read weaklyCoupled3DofFsiFsiDict header
            IOobject weaklyCoupled3DofFsiFsiHeader
            (
                "weaklyCoupled3DofFsiFsiDict",
                yDispl.mesh().time().timeName(),
                "uniform",
                yDispl.mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );
            
            if (weaklyCoupled3DofFsiFsiHeader.headerOk())
            {
                weaklyCoupled3DofFsiFsiDictPtr.reset
                (
                    new IOdictionary
                    (
                        weaklyCoupled3DofFsiFsiHeader
                    )
                );
                
                Info << "Old state restored" << endl;
                weaklyCoupled3DofFsiFsiDictPtr().lookup("YOld") >> Y_;
                Yold_ = Y_;
            }
            
            setDisplacements(yDispl);
    }
    
    return true;
}
void Foam::functionObjects::weaklyCoupled3DofFsi::setDisplacements(volVectorField& yDispl)
{
    if (coordSys_.empty())
    {
        return;
    }

    if (Pstream::parRun())
    {
        Pstream::scatter<vector>(Y_.first());
    }
    
    vector YPatch (coordSys_().globalVector(Y_.first())); //displacements are relative to initial position
    
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchId = iter.key();
        forAll(yDispl.boundaryField()[patchId], faceI)
        {
            yDispl.boundaryFieldRef()[patchId][faceI] = YPatch;
        }
    }
}

bool Foam::functionObjects::weaklyCoupled3DofFsi::execute()
{
    return true;
}

bool Foam::functionObjects::weaklyCoupled3DofFsi::write()
{
    if (!forces::write())
    {
        return false;
    }
    
    volVectorField& yDispl = 
        const_cast<volVectorField&>
        (
            obr_.lookupObject<volVectorField>("cellDisplacement")
        );
    
    if (Pstream::master())
    {
        scalar dt = yDispl.mesh().time().deltaT().value();
        scalar ct = yDispl.mesh().time().value();
        
        vector force = forceEff();
        vector yForce = coordSys_().localVector(force); //convert force to local coord system
        
        //Runge-Kutta 2-nd order method
        Pair<vector> Ymid;
        
        forAll(Ymid.first(), iCmpt)
        {
            Ymid.first()[iCmpt] = Yold_.first()[iCmpt] + 0.5*dt*Yold_.second()[iCmpt];
            Ymid.second()[iCmpt]= Yold_.second()[iCmpt] + 
                0.5*dt*
                (
                    - C_[iCmpt]*Yold_.second()[iCmpt]
                    - K_[iCmpt]*Yold_.first()[iCmpt]
                    + R_*yForce[iCmpt]
                ) / M_[iCmpt];
        
            Y_.first()[iCmpt] = Yold_.first()[iCmpt] + dt*Ymid.second()[iCmpt];
            Y_.second()[iCmpt]= Yold_.second()[iCmpt] + 
                dt*
                (
                    - C_[iCmpt]*Ymid.second()[iCmpt] 
                    - K_[iCmpt]*Ymid.first()[iCmpt] 
                    + R_*yForce[iCmpt]
                ) / M_[iCmpt];
        
            if (mag(Y_.first()[iCmpt]) >= Ymax_[iCmpt])
            {
                Y_.first()[iCmpt] = sign(Y_.first()[iCmpt])*Ymax_[iCmpt];
                Y_.second()[iCmpt] = (Y_.first()[iCmpt] - Yold_.first()[iCmpt]) / dt;
            }
        }

        Yold_ = Y_;
        
        Log << "yForce = " << yForce << endl;
        Log << "Y= " << Y_.first() << endl;
        Log << "Vy= " << Y_.second() << endl;
        
        file(2)<< ct;
        forAll(Y_.first(), iCmpt)
        {
             file(2)<< ";" << Y_.first()[iCmpt];
        }
        forAll(Y_.second(), iCmpt)
        {
            file(2)<< ";" << Y_.second()[iCmpt];
        }
        forAll(yForce, iCmpt)
        {
            file(2)<< ";" << yForce[iCmpt];
        }
        file(2)<< endl;
    }
    
    Pstream::scatter<Pair<vector> >(Yold_);
    
    //write data to file if time is equal to output time
    if (yDispl.mesh().time().outputTime())
    {
        IOdictionary weaklyCoupled3DofFsiFsiDict
        (
            IOobject
            (
                "weaklyCoupled3DofFsiFsiDict",
                yDispl.mesh().time().timeName(),
                "uniform",
                yDispl.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        
        weaklyCoupled3DofFsiFsiDict.set<Pair<vector> >
        (
            "YOld",
            Yold_
        );
        
        weaklyCoupled3DofFsiFsiDict.regIOobject::write();
    }
    
    setDisplacements(yDispl);
    
    return true;
}



// ************************************************************************* //
