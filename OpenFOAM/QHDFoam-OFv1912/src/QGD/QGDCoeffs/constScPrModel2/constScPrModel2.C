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
    This file is part of QGDsolver library, based on OpenFOAM+.

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

#include "constScPrModel2.H"
#include "QGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(constScPrModel2,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        constScPrModel2,
        dictionary
    );
}
}

Foam::qgd::
constScPrModel2::constScPrModel2
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict)
{
    scalar ScQGD = 1.0, PrQGD = 1.0;

    dict.lookup("ScQGD") >> ScQGD;
    dict.lookup("PrQGD") >> PrQGD;

    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;

    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
}

Foam::qgd::
constScPrModel2::~constScPrModel2()
{
}

void Foam::qgd::
constScPrModel2::correct(const Foam::QGDThermo& qgdThermo)
{
    const volScalarField& cSound = qgdThermo.c();
    const volScalarField& p      = qgdThermo.p();
    const volScalarField& mu     = qgdThermo.mu();

    this->tauQGDf_= linearInterpolate(this->aQGD_ / cSound) * hQGDf_;
    this->tauQGD_ = this->aQGD_ * this->hQGD_  / cSound;

    forAll(p.primitiveField(), celli)
    {
        muQGD_.primitiveFieldRef()[celli] =
            p.primitiveField()[celli] *
            ScQGD_.primitiveField()[celli] *
            tauQGD_.primitiveField()[celli];

        alphauQGD_.primitiveFieldRef()[celli] = muQGD_.primitiveField()[celli] /
            PrQGD_.primitiveField()[celli];
    }

    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            muQGD_.boundaryFieldRef()[patchi][facei] =
                p.boundaryField()[patchi][facei] *
                ScQGD_.boundaryField()[patchi][facei] *
                tauQGD_.boundaryField()[patchi][facei];

            alphauQGD_.boundaryFieldRef()[patchi][facei] =
                muQGD_.boundaryFieldRef()[patchi][facei] /
                PrQGD_.boundaryField()[patchi][facei];
        }
    }

    //update tau with molecular viscosity
    this->tauQGD_ += mu / (p * this->ScQGD_);
}

//
//END-OF-FILE
//
