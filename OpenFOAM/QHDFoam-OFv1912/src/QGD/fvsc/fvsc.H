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
Group grpfvsc
    This group contains common part of QGD solvers.
Class
    Foam::fvsc
Description 
    Methods calculating of differential operators
Source files
    fvsc.C   
\*---------------------------------------------------------------------------*/
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "tmp.H"

#ifndef fvsc_H
#define fvsc_H

namespace Foam
{

namespace fvsc
{

    //returns gradient of volScalarField at the surface mesh
    tmp<surfaceVectorField> grad(const volScalarField& vF);
    
    //-
    tmp<surfaceVectorField> grad(const tmp<volScalarField>& vF);

    //returns gradient of volVectorField at the surface mesh
    tmp<surfaceTensorField> grad(const volVectorField& vF);
    
    //-
    tmp<surfaceTensorField> grad(const tmp<volVectorField>& vF);

    //returns divergence of volVectorField at the surface mesh
    tmp<surfaceScalarField> div(const volVectorField& vF);
    
    //-
    tmp<surfaceScalarField> div(const tmp<volVectorField>& vF);

    //returns divergence of volTensorField at the surface mesh
    tmp<surfaceVectorField> div(const volTensorField& vF);
    
    //-
    tmp<surfaceVectorField> div(const tmp<volTensorField>& vF);

} //namespace fvsc

} //namespace Foam

#endif

//END-OF-FILE

