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

#include "simpleFoamLibrary.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::simpleFoamLibrary>
Foam::simpleFoamLibrary::New
(
    const fvMesh& mesh
)
{
    IOdictionary simpleFoamLibraryDict
    (
        IOobject
        (
            "simpleFoamLibraryDict",
            mesh.time().constant(),
            mesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word simpleFoamLibraryTypeName
    (
        simpleFoamLibraryDict.lookup("simpleFoamLibrary")
    );

    Info<< "Selecting simpleFoamLibrary "
        << simpleFoamLibraryTypeName << endl;

    componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_
            ->find(simpleFoamLibraryTypeName);

    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "simpleFoamLibrary::New"
        )   << "Unknown simpleFoamLibrary type "
            << simpleFoamLibraryTypeName << endl << endl
            << "Valid  simpleFoamLibrarys are : " << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<simpleFoamLibrary>(cstrIter()(mesh));
}


// ************************************************************************* //
