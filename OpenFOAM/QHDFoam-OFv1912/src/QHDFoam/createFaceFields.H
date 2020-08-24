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
    createFaceFields
Description
    Create Face Fields for calculation
SourceFile
    QHDFoam.C
\*---------------------------------------------------------------------------*/
     /*
     *
     * Linear interpolation of fields from volumes to face centers
     *
     */
    // Density
    surfaceScalarField rhof
    (
        "rhof",
        linearInterpolate(rho)
    );
    
    // Velocity
    surfaceVectorField Uf
    (
        "Uf",
        linearInterpolate(U)
    );

    surfaceScalarField pf
    (
        "pf",
	linearInterpolate(p)
    );
    
    //
    surfaceVectorField Wf
    (
        "Wf",
        linearInterpolate(W)
    );
    
    surfaceVectorField BdFrcf
    (
        "BdFrcf",
        linearInterpolate(BdFrc)
    );
    
    surfaceScalarField Tf
    (
        "Tf",
	linearInterpolate(T)
    );
    
    // Heat capacity at constant pressure
    surfaceScalarField Cpf
    (
        "Cpf",
        linearInterpolate(thermo.Cp())
    );
    
    surfaceScalarField muf
    (
        "muf",
        linearInterpolate(thermo.mu()*0.0)
    );
    
    /*
     *
     * QGD coefficients
     *
     */
    // OpenFOAM delta coeffs can be negative for coupled patches
    surfaceScalarField taubyrhof
    (
        "taubyrhof",
        tauQGDf/rhof
    );
    
    surfaceScalarField alphaf
    (
        "alphaf",
        linearInterpolate(thermo.alpha())
    );


