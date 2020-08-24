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
    This file is part of QGDsolver, based on OpenFOAM library.

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
    grpPsiQGDThermo
\*---------------------------------------------------------------------------*/



#include "psiQGDThermo.H"
#include "surfaceFields.H"
#include "typeInfo.H"
#include "coupledPolyPatch.H"
#include "coupledFvPatch.H"
#include "emptyFvPatch.H"
#include "wedgeFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "symmetryFvPatch.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiQGDThermo, 0);
    defineRunTimeSelectionTable(psiQGDThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiQGDThermo::psiQGDThermo(const fvMesh& mesh, const word& phaseName)
:
    psiThermo(mesh, phaseName),
    QGDThermo(mesh, *this),
    c_
    (
        "thermo:c",
        qgdCoeffs().hQGD()/mesh.time().deltaT()
    ),
    gamma_
    (
        "thermo:gamma",
        c_ / c_
    )
{
    this->read();
}

Foam::psiQGDThermo::psiQGDThermo(const fvMesh& mesh, const word& phaseName, const word& dictName)
:
    psiThermo(mesh, phaseName, dictName),
    QGDThermo(mesh, *this),
    c_
    (
        "thermo:c",
        qgdCoeffs().hQGD()/mesh.time().deltaT()
    ),
    gamma_
    (
        "thermo:gamma",
        c_ / c_
    )
{
    this->read();
}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiQGDThermo> Foam::psiQGDThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<psiQGDThermo>(mesh, phaseName);
}

Foam::autoPtr<Foam::psiQGDThermo> Foam::psiQGDThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    return basicThermo::New<psiQGDThermo>(mesh, phaseName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiQGDThermo::~psiQGDThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



bool Foam::psiQGDThermo::read()
{
    if (!basicThermo::read())
    {
        return false;
    }

    if (!QGDThermo::read())
    {
      return false;
    }

    return true;
}

const Foam::volScalarField& Foam::psiQGDThermo::c() const
{
    return this->c_;
}

const Foam::volScalarField& Foam::psiQGDThermo::p() const
{
  return psiThermo::p();
}

Foam::volScalarField& Foam::psiQGDThermo::p()
{
  return psiThermo::p();
}

Foam::tmp<Foam::volScalarField> Foam::psiQGDThermo::rho() const
{
  return psiThermo::rho();
}

Foam::tmp<Foam::volScalarField> Foam::psiQGDThermo::mu() const
{
  return psiThermo::mu();
}

// ************************************************************************* //
