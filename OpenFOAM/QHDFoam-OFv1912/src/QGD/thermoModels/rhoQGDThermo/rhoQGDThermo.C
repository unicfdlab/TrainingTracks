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
    grpRhoQGDThermo

\*---------------------------------------------------------------------------*/

#include "rhoQGDThermo.H"
#include "QGDCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rhoQGDThermo, 0);
    defineRunTimeSelectionTable(rhoQGDThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoQGDThermo::rhoQGDThermo(const fvMesh& mesh, const word& phaseName)
:
    rhoThermo(mesh, phaseName),
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

Foam::rhoQGDThermo::rhoQGDThermo(const fvMesh& mesh, const word& phaseName, const word& dictName)
:
    rhoThermo(mesh, phaseName,dictName),
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

Foam::autoPtr<Foam::rhoQGDThermo> Foam::rhoQGDThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<rhoQGDThermo>(mesh, phaseName);
}

Foam::autoPtr<Foam::rhoQGDThermo> Foam::rhoQGDThermo::New
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictName
)
{
    return basicThermo::New<rhoQGDThermo>(mesh, phaseName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoQGDThermo::~rhoQGDThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



bool Foam::rhoQGDThermo::read()
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

const Foam::volScalarField& Foam::rhoQGDThermo::c() const
{
    return this->c_;
}

const Foam::volScalarField& Foam::rhoQGDThermo::p() const
{
  return rhoThermo::p();
}

Foam::volScalarField& Foam::rhoQGDThermo::p()
{
  return rhoThermo::p();
}

Foam::tmp<Foam::volScalarField> Foam::rhoQGDThermo::rho() const
{
  return rhoThermo::rho();
}

Foam::tmp<Foam::volScalarField> Foam::rhoQGDThermo::mu() const
{
  return rhoThermo::mu();
}

// ************************************************************************* //
