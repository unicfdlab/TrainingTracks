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

#include "weaklyCoupledFsi.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "wordReList.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "dynamicFvMesh.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(weaklyCoupledFsi, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::weaklyCoupledFsi::weaklyCoupledFsi
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    forces
    (
	name,
	obr,
	dict,
	loadFromFiles,
	readFields
    ),
    M_(0.0),
    C_(0.0),
    K_(0.0),
    R_(0.0),
    Ymax_(0.0),
    append_(false),
    Y_ (0.0, 0.0),
    Yold_(0.0, 0.0),
    out_(NULL)
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

        out_.reset
        (
            new OFstream
            (
                dict.lookup("results")
            )
        );
    
        if (append_ && oldFileLines.size())
        {
            forAll(oldFileLines, iLine)
            {
                out_() << oldFileLines[iLine] << endl;
            }
        }
        else
        {
            out_() << "Time;Y;Vy;Fy" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::weaklyCoupledFsi::~weaklyCoupledFsi()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::weaklyCoupledFsi::read(const dictionary& dict)
{
    forces::read(dict);
    
    dict.lookup("M") >> M_;
    
    dict.lookup("C") >> C_;
    
    dict.lookup("K") >> K_;
    
    dict.lookup("R") >> R_;
    
    dict.lookup("Ymax") >> Ymax_;
    
    dict.lookup("append") >> append_;

    Info << "Reading old state" << endl;
    autoPtr<IOdictionary> weaklyCoupledFsiDictPtr;
    //try to read weaklyCoupledFsi object properties
    {
        volVectorField& yDispl = 
            const_cast<volVectorField&>
            (
                obr_.lookupObject<volVectorField>("cellDisplacement")
            );
    
        //read weaklyCoupledFsiDict header
        IOobject weaklyCoupledFsiHeader
        (
            "weaklyCoupledFsiDict",
            yDispl.mesh().time().timeName(),
            yDispl.mesh(),
            IOobject::MUST_READ
        );
        
        if (weaklyCoupledFsiHeader.headerOk())
        {
            weaklyCoupledFsiDictPtr.reset
            (
                new IOdictionary
                (
                    weaklyCoupledFsiHeader
                )
            );
        
            weaklyCoupledFsiDictPtr().lookup("YOld") >> Y_;
            Yold_ = Y_;
        }
        
        setDisplacements(yDispl);
    }

}

void Foam::weaklyCoupledFsi::setDisplacements(volVectorField& yDispl)
{
    if (Pstream::parRun())
    {
        Pstream::scatter<scalar>(Y_.first());
    }

    vector YPatch (0.0, Y_.first(), 0.0);
    
    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchId = iter.key();
        forAll(yDispl.boundaryField()[patchId], faceI)
        {
            yDispl.boundaryField()[patchId][faceI] = YPatch;
        }
    }
}

void Foam::weaklyCoupledFsi::write()
{

    if (!active_)
    {
        return;
    }

    forces::write();

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
        scalar yForce = force.y();
        
        Pair<scalar> Ymid;
        
        Ymid.first() = Yold_.first() + 0.5*dt*Yold_.second();
        Ymid.second()= Yold_.second() + 0.5*dt*(-C_*Yold_.second() - K_*Yold_.first() + R_*yForce) / M_;
        
        Y_.first() = Yold_.first() + dt*Ymid.second();
        Y_.second()= Yold_.second() + dt*(-C_*Ymid.second() - K_*Ymid.first() + R_*yForce) / M_;
        
        if (mag(Y_.first()) >= Ymax_)
        {
            Y_.first() = sign(Y_.first())*Ymax_;
            Y_.second() = (Y_.first() - Yold_.first()) / dt;
        }
        
        Yold_ = Y_;
        
        if (log_)
        {
            Info << "yForce = " << yForce << endl;
            Info << "Y= " << Y_.first() << endl;
            Info << "Vy= " << Y_.second() << endl;
        }
        
        out_() << ct << ";" << Y_.first() << ";" << Y_.second() << ";" << yForce << endl;
        
        //write data to file if time is equal to output time
        if (yDispl.mesh().time().outputTime())
        {
            IOdictionary weaklyCoupledFsiDict
            (
                IOobject
                (
                    "weaklyCoupledFsiDict",
                    yDispl.mesh().time().timeName(),
                    "uniform",
                    yDispl.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
        
            weaklyCoupledFsiDict.set<Pair<scalar> >
            (
                "YOld",
                Yold_
            );

            weaklyCoupledFsiDict.regIOobject::write();
        }
    }
    
    setDisplacements(yDispl);
}


// ************************************************************************* //
