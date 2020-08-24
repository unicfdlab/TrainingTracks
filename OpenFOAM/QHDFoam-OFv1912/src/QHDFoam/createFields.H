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
Global
    createFields
Description
    Creating fields for calculations
SourceFile
    QHDFoam.C
\*---------------------------------------------------------------------------*/

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoQGDThermo> pThermo
(
    rhoQGDThermo::New(mesh)
);
rhoQGDThermo& thermo = pThermo();
thermo.correct();

volScalarField& e = thermo.he();

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const surfaceScalarField& hQGDf = thermo.hQGDf();
const surfaceScalarField& tauQGDf = thermo.tauQGDf();

Info << "Thermo corrected" << endl;

autoPtr<compressible::turbulenceModel> turbulence;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volVectorField W
(
    IOobject
    (
        "W",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    U
);




Info<< "\nReading gravitationalProperties" << endl;
  
IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g(gravitationalProperties.lookup("g"));

dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    thermo.subDict("mixture").subDict("transport")
);


surfaceScalarField phiu
(
    "phiu",
    mesh.Sf() & linearInterpolate(U)
);

surfaceScalarField  phiwo
(
    "phiwo",
    mesh.Sf() & linearInterpolate(W)
);

surfaceScalarField phi
(
    "phi",
    mesh.Sf() & (linearInterpolate(U) - linearInterpolate(W))
);

volVectorField BdFrc
(
    "BdFrc",
    T*g*beta
);

Switch implicitDiffusion (thermo.implicitDiffusion());

surfaceScalarField phiRhof
(
    "phiRhof",
    linearInterpolate(rho)*phi
);

Info<< "Creating turbulence model\n" << endl;
turbulence.reset
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phiRhof,
        thermo
    ).ptr()
);

label pRefCell = 0;
scalar pRefValue = 0.0;

setRefCell(p, thermo.subDict("QGD"), pRefCell, pRefValue);

#include "createIcoZeroSources.H"

//
//END-OF-FILE
//

