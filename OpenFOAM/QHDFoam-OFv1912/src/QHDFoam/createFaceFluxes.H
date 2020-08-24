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
    createFaceFluxes
Description
    Create Face Fluxes for calculation
SourceFile
    QHDFoam.C
\*---------------------------------------------------------------------------*/
// Field for debugging parallel execution and tau-terms evaluation
/*
volScalarField cellNo
(
    IOobject
    (
        "cellNo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*0.0
);

forAll(cellNo, celli)
{
    cellNo[celli] = mesh.C()[celli].z();
}
forAll(cellNo.boundaryField(), patchi)
{
    forAll(cellNo.boundaryField()[patchi], facei)
    {
        cellNo.boundaryFieldRef()[patchi][facei] = mesh.Cf().boundaryField()[patchi][facei].z();
    }
}
surfaceVectorField cellNoGradf
(
    "cellNoGradf", fvsc::grad(cellNo)
);

cellNo.write();
cellNoGradf.write();
return 0;
*/
//Gradients and divergence
//---------Start---------


surfaceVectorField gradPf 
(
    "gradPf", fvsc::grad(p)
);

surfaceVectorField gradTf
(
    "gradTf", fvsc::grad(T)
);

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceTensorField gradWf
(
    "gradWf",
    fvsc::grad(W)
);

surfaceVectorField gradRhof
(
    "gradRhof",
    fvsc::grad(rho)
);

surfaceScalarField divUf
(
    "divUf",
    tr(gradUf)
);

//---------End---------
//P equation fluxes
//---------Start---------
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
surfaceVectorField phiUf
(
    "phiUf",
    phi * Uf + (mesh.Sf() & (Wf * Wf))
);

surfaceScalarField Pbyrho
(
    "Pbyrho",
    pf/rhof
);

surfaceScalarField phiTf
(
    "phiTf",
    phi*Tf
);

surfaceVectorField phiPi
(
    "phiPi",
    muf/rhof * mesh.Sf() & (gradUf+gradUf.T())
);

//---------End---------

// Fluxes for T balance equation
//---------Start---------
surfaceScalarField phiUmWT
(
    "phiUmWT",
    ((Uf-Wf)*Tf) & mesh.Sf()
);

surfaceScalarField Hif
(
    "Hif",
    alphaf/rhof
);
Info << "Creating Face Fluxes is done"<< endl;

//---------End---------
