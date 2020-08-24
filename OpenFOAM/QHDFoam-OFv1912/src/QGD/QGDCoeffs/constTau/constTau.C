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

#include "constTau.H"
#include "QGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(constTau,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        constTau,
        dictionary
    );
}
}

Foam::qgd::
constTau::constTau
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict),
    tau_(0.0)
{
    scalar ScQGD = 0.0, PrQGD = 1.0;

    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;
    muQGD_.primitiveFieldRef() = 0.0;
    alphauQGD_.primitiveFieldRef() = 0.0;

    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
    muQGD_.boundaryFieldRef() = 0.0;
    alphauQGD_.boundaryFieldRef() = 0.0;
    
    dict.lookup("Tau") >> tau_;

    this->tauQGD_ = dimensionedScalar("tauQGD",dimTime,tau_);
    this->tauQGDf_= linearInterpolate(this->tauQGD_);
}

Foam::qgd::
constTau::~constTau()
{
}

void Foam::qgd::
constTau::correct(const Foam::QGDThermo& qgdThermo)
{
}

//
//END-OF-FILE
//
