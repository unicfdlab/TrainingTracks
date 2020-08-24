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
Group grpreduced
    This group contains common part of QGD solvers.
Class
    Foam::fvsc::reduced::reducedFaceNormalStencil
Description 
    Methods calculating of differential operators without tangential direvatives  
\*---------------------------------------------------------------------------*/



#include "fvc.H"
#include "reducedFaceNormalStencil.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fvsc
{
    defineTypeNameAndDebug(reduced,0);
    addToRunTimeSelectionTable
    (
        fvscStencil,
        reduced,
        components
    );
}
}

// constructors
Foam::fvsc::reduced::reduced(const IOobject& io)
:
    fvscStencil(io)
{
}

Foam::fvsc::reduced::~reduced()
{
}

//- Calculate gradient of volume scalar function on the faces
//
// \param iF         Internal scalar field.
//                   Allowable values: constant reference to the volScalarField.
//
// \return           Gradient of iF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::reduced::Grad(const volScalarField& vF)
{
    tmp<surfaceVectorField> tgradIF(nf_ * fvc::snGrad(vF));
    
    return tgradIF;
};

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceTensorField> Foam::fvsc::reduced::Grad(const volVectorField& iVF)
{

    tmp<surfaceTensorField> tgradIVF(nf_ * fvc::snGrad(iVF));

    return tgradIVF;
};

Foam::tmp<Foam::surfaceScalarField> Foam::fvsc::reduced::Div(const volVectorField& iVF)
{
    tmp<surfaceScalarField> tdivIVF(nf_ & fvc::snGrad(iVF));
    
    return tdivIVF;
};

//- Calculate divergence of volume tensor field on the faces.
//
// \param iTF        Internal tensor field.
//                   Allowable values: constant reference to the volTensorField.
//
// \return           Divergence of iTF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::reduced::Div(const volTensorField& iTF)
{
    tmp<surfaceVectorField> tdivITF(nf_ & fvc::snGrad(iTF));
    
    return tdivITF;
}


//
//END-OF-FILE
//


