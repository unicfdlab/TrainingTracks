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
Group 
    grpGaussVolPoint
Class
    Foam::fvsc::GaussVolPoint::GaussVolPointStencil    
Description
    Methods calculation of differential operators without tangential derivatives
    
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "GaussVolPointStencil.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"
#include "processorFvPatchField.H"

namespace Foam
{
namespace fvsc
{
    defineTypeNameAndDebug(GaussVolPoint,0);
    addToRunTimeSelectionTable
    (
        fvscStencil,
        GaussVolPoint,
        components
    );
}
}

// constructors
Foam::fvsc::GaussVolPoint::GaussVolPoint(const IOobject& io)
:
    fvscStencil(io),
    GaussVolPointBase(refCast<const fvMesh>(io.db()))
{
}

Foam::fvsc::GaussVolPoint::~GaussVolPoint()
{
}

//- Calculate gradient of volume scalar function on the faces
//
// \param iF         Internal scalar field.
//                   Allowable values: constant reference to the volScalarField.
//
// \return           Gradient of iF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::GaussVolPoint::Grad(const volScalarField& vF)
{
    const_cast<volScalarField&>(vF).correctBoundaryConditions();
    tmp<surfaceVectorField> tgradIF(vector::zero * fvc::snGrad(vF));
    //tmp<surfaceVectorField> tgradIF(nf_ * fvc::snGrad(vF));
    surfaceVectorField& gradIf = tgradIF.ref();
    
    this->faceGrad(vF, gradIf);
    
    return tgradIF;
};

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceTensorField> Foam::fvsc::GaussVolPoint::Grad(const volVectorField& iVF)
{
    const_cast<volVectorField&>(iVF).correctBoundaryConditions();
    tmp<surfaceTensorField> tgradIVF(tensor::zero * fvc::snGrad(iVF.component(0)));
    //tmp<surfaceTensorField> tgradIVF(nf_ * fvc::snGrad(iVF));
    surfaceTensorField& gradIVF = tgradIVF.ref();
    
    this->faceGrad(iVF, gradIVF);
    
    return tgradIVF;
};

Foam::tmp<Foam::surfaceScalarField> Foam::fvsc::GaussVolPoint::Div(const volVectorField& iVF)
{
    const_cast<volVectorField&>(iVF).correctBoundaryConditions();
    tmp<surfaceScalarField> tdivIVF(0.0 * fvc::snGrad(iVF.component(0)));
    //tmp<surfaceScalarField> tdivIVF(nf_ & fvc::snGrad(iVF));
    surfaceScalarField& divIVF = tdivIVF.ref();
    
    this->faceDiv(iVF, divIVF);
    
    return tdivIVF;
};

//- Calculate divergence of volume tensor field on the faces.
//
// \param iTF        Internal tensor field.
//                   Allowable values: constant reference to the volTensorField.
//
// \return           Divergence of iTF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::GaussVolPoint::Div(const volTensorField& iTF)
{
    const_cast<volTensorField&>(iTF).correctBoundaryConditions();
    tmp<surfaceVectorField> tdivITF(vector::zero*fvc::snGrad(iTF.component(0)));
    //tmp<surfaceVectorField> tdivITF(nf_*fvc::snGrad(iTF.component(0)));
    surfaceVectorField& divITF = tdivITF.ref();

    this->faceDiv(iTF, divITF);

    return tdivITF;
}


//
//END-OF-FILE
//


