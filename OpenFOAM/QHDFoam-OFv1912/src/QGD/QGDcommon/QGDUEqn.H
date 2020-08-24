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
    grpQGDcommon

Description
    Solution of momentum equation for QGD solver

\*----------------------------------------------------------------------------*/
    solve
    (
        fvm::ddt(rhoU)
        +
        fvc::div(phiJmU)
        +
        fvc::div(phiP)
        -
        fvc::div(phiPi)
    );
    
    // Correct velocity
    U.ref() =
        rhoU()
       /rho();
    U.correctBoundaryConditions();
    
   // Solve diffusive QGD & NS part
    if (implicitDiffusion)
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U) - fvc::ddt(rho,U)
            - fvm::laplacian(muf, U)
            - fvc::div(phiTauMC)
            ==
            rhoUSu
        );
        
        solve
        (
            UEqn
        );
        
        rhoU = rho*U;
        
        sigmaDotUPtr() = (muf*linearInterpolate(fvc::grad(U)) + tauMCPtr()) & Uf;
        
        phiSigmaDotU = mesh.Sf() & sigmaDotUPtr(); //or eqn.flux()?
    }
    else
    {
        // Correct velocity with distributed sources
        solve
        (
            fvm::ddt(rho,U)
            -
            fvc::ddt(rhoU)
            ==
            rhoUSu
        );
    }
    rhoU.boundaryFieldRef() == rho.boundaryField()*
    U.boundaryField();
    
//
//END-OF-FILE
//


