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

Class
    Foam::altSimpleFoamLibrary

Description

SourceFiles
    altSimpleFoamLibrary.C


\*---------------------------------------------------------------------------*/

#ifndef altSimpleFoamLibrary_H
#define altSimpleFoamLibrary_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "autoPtr.H"
#include "simpleFoamLibrary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class altSimpleFoamLibrary Declaration
\*---------------------------------------------------------------------------*/

class altSimpleFoamLibrary
:
    public simpleFoamLibrary
{

protected:

    // Protected data

        //- Disallow copy construct
        altSimpleFoamLibrary(const altSimpleFoamLibrary&);

        //- Disallow default bitwise assignment
        void operator=(const altSimpleFoamLibrary&);

public:

    //- Runtime type information
    TypeName("altSimpleFoamLibrary");

    // Constructors

        //- Construct from components
        altSimpleFoamLibrary
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~altSimpleFoamLibrary()
    {}


    // Member Functions
    virtual void execute();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
