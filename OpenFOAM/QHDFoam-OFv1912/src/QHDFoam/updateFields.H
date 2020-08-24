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
    updateFields
Description
    Update fields during calculations
SourceFile
    QHDFoam.C
\*---------------------------------------------------------------------------*/
/*
 *
 * Gradients at faces
 *
 */
gradUf = fvsc::grad(U);

gradWf = fvsc::grad(W);

gradTf = fvsc::grad(T);

/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
muf    = qgdInterpolate(turbulence->muEff());
alphaf = qgdInterpolate(thermo.alpha());

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = qgdInterpolate(rho);

// Velocity
Uf   = qgdInterpolate(U);
    
// T
Tf   = qgdInterpolate(T);

//
BdFrc  = beta*T*g;
BdFrcf = qgdInterpolate(BdFrc);

//Hi
Hif = alphaf/rhof;

// Heat capacity at constant pressure
Cpf = qgdInterpolate(thermo.Cp());




