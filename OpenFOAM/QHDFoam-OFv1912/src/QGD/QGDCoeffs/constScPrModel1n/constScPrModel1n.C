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

#include "constScPrModel1n.H"
#include "QGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(constScPrModel1n,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        constScPrModel1n,
        dictionary
    );
}
}

Foam::qgd::
constScPrModel1n::constScPrModel1n
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict)
{
    scalar PrQGD = 1.0;
    if (dict.found("PrQGD"))
    {
        dict.lookup("PrQGD") >> PrQGD;
    }
    PrQGD_.primitiveFieldRef() = PrQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;

    IOobject ScHeader
    (
        "ScQGD",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );

    if (ScHeader.typeHeaderOk<volScalarField>())
    {
        //do nothing, ScQGD field is present
        ScQGD_.writeOpt() = IOobject::AUTO_WRITE;
    }
    else
    {
        scalar ScQGD = 1.0;
        if (dict.found("ScQGD"))
        {
            dict.lookup("ScQGD") >> ScQGD;
        }
        ScQGD_.primitiveFieldRef() = ScQGD;
        ScQGD_.boundaryFieldRef() = ScQGD;
    }
}

Foam::qgd::
constScPrModel1n::~constScPrModel1n()
{
}

void Foam::qgd::
constScPrModel1n::correct(const Foam::QGDThermo& qgdThermo)
{
    const volScalarField& cSound = qgdThermo.c();
    const volScalarField& p      = qgdThermo.p();
    if (!p.mesh().thisDb().foundObject<volVectorField>("U"))
    {
        this->tauQGDf_= linearInterpolate(this->aQGD_) * hQGDf_ /linearInterpolate(cSound);
        this->tauQGD_ = this->aQGD_ * this->hQGD_  / cSound;
    }
    else
    {
        const surfaceScalarField magUn
        (
            mag
            (
                linearInterpolate
                (
                   p.mesh().thisDb().lookupObject<volVectorField>("U")
                ) & (p.mesh().Sf()/p.mesh().magSf())
            )
        );
        const volScalarField magU
        (
            mag
            (
               p.mesh().thisDb().lookupObject<volVectorField>("U")
            )
        );
        this->tauQGD_ = this->aQGD_ * this->hQGD_  / (magU + cSound);
        //this->tauQGDf_= linearInterpolate(this->aQGD_) * hQGDf_ /(magUn + linearInterpolate(cSound));
        this->tauQGDf_= linearInterpolate(this->tauQGD_);
    }
    
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
}

//
//END-OF-FILE
//
