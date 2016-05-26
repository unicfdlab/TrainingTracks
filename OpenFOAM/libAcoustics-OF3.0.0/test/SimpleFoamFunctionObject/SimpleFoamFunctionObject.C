/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "SimpleFoamFunctionObject.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(SimpleFoamFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, SimpleFoamFunctionObject, dictionary);
}

// * * * * * * * * * * * * * * * Private Members * * * * * * * * * * * * * * //

void Foam::SimpleFoamFunctionObject::readDict()
{
    dict_.lookup("helloAddress") >> helloAddress_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SimpleFoamFunctionObject::SimpleFoamFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    dict_(dict)
{
    readDict();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::SimpleFoamFunctionObject::start()
{
    readDict();
    
    Info << "SimpleFoamFunctionObject::start() method tells: Hello, " << helloAddress_ << ", current problem time is: " << time_.value() << endl;

    return true;
}

bool Foam::SimpleFoamFunctionObject::execute(bool)
{
    readDict();
    
    Info << "SimpleFoamFunctionObject::execute() method tells: Hello, " << helloAddress_ << ", current problem time is: " << time_.value() << endl;

    return true;
}

bool Foam::SimpleFoamFunctionObject::end()
{
    readDict();
    
    Info << "SimpleFoamFunctionObject::end() method tells: Hello, " << helloAddress_ << ", current problem time is: " << time_.value() << endl;

    return true;
}

bool Foam::SimpleFoamFunctionObject::read
(
    const dictionary& dict
)
{
    if (dict != dict_)
    {
        dict_ = dict;

        return start();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
