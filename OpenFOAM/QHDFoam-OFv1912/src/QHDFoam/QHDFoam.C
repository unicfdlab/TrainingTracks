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
Application
    QHDFoam
Description
    Solver for unsteady 3D turbulent flow of incompressible viscous fluid 
    governed by quasi-hydrodynamic dynamic (QHD) equations.

    QHD system of equations has been developed by scientific group from
    Keldysh Institute of Applied Mathematics,
    see http://elizarova.imamod.ru/selection-of-papers.html
    A comprehensive description of QGD equations and their applications
    can be found here:
    \verbatim
    Elizarova, T.G.
    "Quasi-Gas Dynamic equations"
    Springer, 2009
    \endverbatim
    A brief of theory on QGD and QHD system of equations:
    \verbatim
    Elizarova, T.G. and Sheretov, Y.V.
    "Theoretical and numerical analysis of quasi-gasdynamic and quasi-hydrodynamic
    equations"
    J. Computational Mathematics and Mathematical Physics, vol. 41, no. 2, pp 219-234,
    2001
    \endverbatim

    Developed by UniCFD group (www.unicfd.ru) of ISP RAS (www.ispras.ru).
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "QHD.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFaceFields.H"
    #include "createFaceFluxes.H"
    #include "createTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        /*
         *
         * Update fields
         *
         */
        #include "updateFields.H"

        /*
         *
         * Update fluxes
         *
         */
       #include "updateFluxes.H"

        /*
         *
         * Update time step
         *
         */
        #include "readTimeControls.H"
        #include "QHDCourantNo.H"
        #include "setDeltaT-QGDQHD.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Store old time values
        U.oldTime();
        T.oldTime();
        turbulence->correct();
        
        #include "QHDpEqn.H"
        
        #include "QHDUEqn.H"
        
        #include "QHDTEqn.H"
        
        if (p.needReference())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pRefValue - getRefCellValue(p, pRefCell)
            );
        }
        
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
    }



    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
