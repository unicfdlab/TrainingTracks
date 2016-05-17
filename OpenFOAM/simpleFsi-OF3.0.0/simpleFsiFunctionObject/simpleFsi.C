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

#include "simpleFsi.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "wordReList.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleFsi, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleFsi::simpleFsi
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
    Y_ (0.0, 0.0),
    Yold_(0.0, 0.0),
    yD_(NULL)
{
    this->read(dict);
    yD_.reset
    (
        new OFstream
        (
            dict.lookup("results")
        )
    );
    
    yD_() << "Time;Y;Vy;Fy" << endl;
}


Foam::simpleFsi::simpleFsi
(
    const word& name,
    const objectRegistry& obr,
    const labelHashSet& patchSet,
    const word& pName,
    const word& UName,
    const word& rhoName,
    const scalar rhoInf,
    const scalar pRef,
    const coordinateSystem& coordSys
)
:
    forces
    (
	name,
	obr,
	patchSet,
	pName,
	UName,
	rhoName,
	rhoInf,
	pRef,
	coordSys
    ),
    M_(0.0),
    C_(0.0),
    K_(0.0),
    R_(0.0),
    Ymax_(0.0),
    Y_ (0.0, 0.0),
    Yold_(0.0, 0.0),
    yD_(NULL)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleFsi::~simpleFsi()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::simpleFsi::read(const dictionary& dict)
{
    forces::read(dict);
    
    dict.lookup("M") >> M_;
    
    dict.lookup("C") >> C_;
    
    dict.lookup("K") >> K_;
    
    dict.lookup("R") >> R_;
    
    dict.lookup("Ymax") >> Ymax_;
}

void Foam::simpleFsi::write()
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
	
	vector force = sum(force_[0] + force_[1]);
	
	scalar yForce = force.y();
	
	Pair<scalar> Ymid;
	
	Ymid.first() = Yold_.first() + 0.5*dt*Yold_.second();
	Ymid.second()= Yold_.second() + 0.5*dt/ M_ * (-C_*Yold_.second() - K_*Yold_.first() + R_*yForce);
	
	Y_.first() = Yold_.first() + dt*Ymid.second();
	Y_.second()= Yold_.second() + dt/ M_ * (-C_*Ymid.second() - K_*Ymid.first() + R_*yForce);
	
	if (mag(Y_.first()) >= Ymax_)
	{
	    Y_.first() = sign(Y_.first())*Ymax_;
	    Y_.second() = (Y_.first() - Yold_.first()) / dt;
	}
	
	Yold_ = Y_;
	
	Info << "yForce = " << yForce << endl;
	Info << "Y= " << Y_.first() << endl;
	Info << "Vy= " << Y_.second() << endl;
	
	yD_() << ct << ";" << Y_.first() << ";" << Y_.second() << ";" << yForce << endl;
	
	if (Pstream::parRun())
	{
	    for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
	    {
		OPstream::write(Pstream::scheduled, jSlave, reinterpret_cast<char*>(&(Y_.first())), sizeof(scalar),
		UPstream::msgType(), UPstream::worldComm);
	    }
	}
    }
    else
    {
	IPstream::read(Pstream::scheduled, Pstream::masterNo(), reinterpret_cast<char*>(&(Y_.first())), sizeof(scalar),
		        UPstream::msgType(), UPstream::worldComm);
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



// ************************************************************************* //
