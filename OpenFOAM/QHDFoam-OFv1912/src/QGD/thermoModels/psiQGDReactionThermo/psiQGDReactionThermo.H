  /*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2016-2019 ISP RAS (www.ispras.ru) UniCFD Group (www.unicfd.ru)
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
    Foam::psiQGDReactionThermo

Group
    grpPsiQGDReactionThermo

Description
    Foam::psiQGDReactionThermo

SourceFiles
    psiQGDReactionThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiQGDReactionThermo_H
#define psiQGDReactionThermo_H

#include "psiQGDThermo.H"
#include "basicSpecieMixture.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class psiQGDReactionThermo Declaration
\*---------------------------------------------------------------------------*/

class psiQGDReactionThermo
:
    public psiQGDThermo
{

public:

    //- Runtime type information
    TypeName("psiQGDReactionThermo");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiQGDReactionThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiQGDReactionThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName)
    );


    // Constructors

        //- Construct from mesh and phase name
        psiQGDReactionThermo
        (
            const fvMesh&,
            const word& phaseName
        );


        //- Construct from mesh,dictionary,phase name with a single temperature
        psiQGDReactionThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<psiQGDReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


        //- Selector based on mesh, phase name and dictionary
        static autoPtr<psiQGDReactionThermo> New
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Destructor
    virtual ~psiQGDReactionThermo();


    // Member functions

        //- Return the composition of the multi-component mixture
        virtual basicSpecieMixture& composition() = 0;

        //- Return the composition of the multi-component mixture
        virtual const basicSpecieMixture& composition() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
