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
    
Group
    grpPsiQGDReactionThermo

\*---------------------------------------------------------------------------*/

#include "psiQGDReactionThermo.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiQGDReactionThermo, 0);
    defineRunTimeSelectionTable(psiQGDReactionThermo, fvMesh);
    defineRunTimeSelectionTable(psiQGDReactionThermo, fvMeshDictPhase);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiQGDReactionThermo::psiQGDReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    psiQGDThermo(mesh, phaseName)
{}


Foam::psiQGDReactionThermo::psiQGDReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
:
    psiQGDThermo(mesh, phaseName, dictName)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiQGDReactionThermo> Foam::psiQGDReactionThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<psiQGDReactionThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::psiQGDReactionThermo> Foam::psiQGDReactionThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    return basicThermo::New<psiQGDReactionThermo>(mesh, phaseName, dictName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiQGDReactionThermo::~psiQGDReactionThermo()
{}


// ************************************************************************* //
