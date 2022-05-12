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
    Foam::simpleFoamLibrary

Description

SourceFiles
    simpleFoamLibrary.C
    newSimpleFoamLibrary.C

\*---------------------------------------------------------------------------*/

#ifndef simpleFoamLibrary_H
#define simpleFoamLibrary_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class simpleFoamLibrary Declaration
\*---------------------------------------------------------------------------*/

class simpleFoamLibrary
{

protected:

    // Protected data
        
        //-
        const fvMesh& mesh_;

        //- Disallow copy construct
        simpleFoamLibrary(const simpleFoamLibrary&);

        //- Disallow default bitwise assignment
        void operator=(const simpleFoamLibrary&);

public:

    //- Runtime type information
    TypeName("simpleFoamLibrary");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            simpleFoamLibrary,
            components,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<simpleFoamLibrary> New
        (
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from components
        simpleFoamLibrary
        (
            const word& type,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~simpleFoamLibrary()
    {}


    // Member Functions
    virtual void execute();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
