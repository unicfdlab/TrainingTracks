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

Group grpQGDcommon
    This group contains common part of QGD solvers.

Description
    Solution of momentum equation for QGD solver.

SourceFiles
    QGDFoam.C

\*----------------------------------------------------------------------------*/
{
    combustion->correct();
    Qdot = combustion->Qdot();
    volScalarField Yt(0.0*Y[0]);
    
    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];
            const surfaceScalarField& phiJmYi = phiJmY[i];
            
            if (implicitDiffusion)
            {
                fvScalarMatrix YEqn
                (
                    fvm::ddt(rho,Yi)
                    +
                    fvc::div(phiJmYi)
                    -
                    fvm::laplacian(muf/ScNumbers[i],Yi)
                    ==
                    combustion->R(Yi)
                    +
                    parcels.SYi(i, Yi)
                );
                
                YEqn.solve();
                
                diffusiveFlux[i] += YEqn.flux();
                diffusiveFlux[inertIndex] -= diffusiveFlux[i];
            }
            else
            {
                solve
                (
                    fvm::ddt(rho,Yi)
                    +
                    fvc::div(phiJmYi)
                    -
                    fvc::laplacian(muf/ScNumbers[i],Yi)
                    ==
                    combustion->R(Yi)
                    +
                    parcels.SYi(i, Yi)
                );
                
                diffusiveFlux[i] += (muf / ScNumbers[i]) * fvc::snGrad(Yi.oldTime()) * mesh.magSf();
                diffusiveFlux[inertIndex] -= diffusiveFlux[i];
            }
            
            Yi.max(0.0);
            Yt += Yi;
        }
        
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
    }
}

//
//END-OF-FILE
//


