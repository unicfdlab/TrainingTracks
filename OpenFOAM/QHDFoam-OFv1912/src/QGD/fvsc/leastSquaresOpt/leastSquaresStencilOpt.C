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
Group grpleastSquaresOpt
    This group contains common part of QGD solvers.
Class
    Foam::fvsc::leastSquaresOpt::leastSquaresStencilOpt
Description 
    Methods for optimal calculating of directional derivative. 
    With parallel realisation.
\*---------------------------------------------------------------------------*/


#include "leastSquaresStencilOpt.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include "wedgeFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "symmetryFvPatch.H"
#include "emptyFvPatch.H"
#include <HashTable.H>

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fvsc
{
    defineTypeNameAndDebug(leastSquaresOpt,0);
    addToRunTimeSelectionTable
    (
        fvscStencil,
        leastSquaresOpt,
        components
    );
}
}

// constructors
Foam::fvsc::leastSquaresOpt::leastSquaresOpt(const IOobject& io)
:
    fvscStencil(io),
    leastSquaresBase(this->mesh_)
{
}


Foam::fvsc::leastSquaresOpt::~leastSquaresOpt()
{
}

Foam::tmp<Foam::surfaceTensorField> Foam::fvsc::leastSquaresOpt::Grad(const volVectorField& iVF)
{
    return Grad(iVF, linearInterpolate(iVF));
}

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceTensorField> Foam::fvsc::leastSquaresOpt::Grad(const volVectorField& iVF, const surfaceVectorField& sVF)
{
    
    tmp<surfaceTensorField> tgradIVF(0*nf_*fvc::snGrad(iVF));
    surfaceScalarField tField = sVF.component(0)*0;
    surfaceTensorField& gradIVF = tgradIVF.ref();

    faceScalarDer(iVF.primitiveField().component(0),sVF.primitiveField().component(0),0,tField);
    gradIVF.primitiveFieldRef().replace(0,tField);

    faceScalarDer(iVF.primitiveField().component(0),sVF.primitiveField().component(0),1,tField);
    gradIVF.primitiveFieldRef().replace(3,tField);

    faceScalarDer(iVF.primitiveField().component(0),sVF.primitiveField().component(0),2,tField);
    gradIVF.primitiveFieldRef().replace(6,tField);

    faceScalarDer(iVF.primitiveField().component(1),sVF.primitiveField().component(1),0,tField);
    gradIVF.primitiveFieldRef().replace(1,tField);

    faceScalarDer(iVF.primitiveField().component(1),sVF.primitiveField().component(1),1,tField);
    gradIVF.primitiveFieldRef().replace(4,tField);

    faceScalarDer(iVF.primitiveField().component(1),sVF.primitiveField().component(1),2,tField);
    gradIVF.primitiveFieldRef().replace(7,tField);

    faceScalarDer(iVF.primitiveField().component(2),sVF.primitiveField().component(2),0,tField);
    gradIVF.primitiveFieldRef().replace(2,tField);

    faceScalarDer(iVF.primitiveField().component(2),sVF.primitiveField().component(2),1,tField);
    gradIVF.primitiveFieldRef().replace(5,tField);

    faceScalarDer(iVF.primitiveField().component(2),sVF.primitiveField().component(2),2,tField);
    gradIVF.primitiveFieldRef().replace(8,tField);


    forAll(mesh_.boundaryMesh(), ipatch)
    {
        bool notConstrain = true;
        const fvPatch& fvp = mesh_.boundary()[ipatch];
        if
        (
            isA<emptyFvPatch>(fvp) ||
            isA<wedgeFvPatch>(fvp) ||
            isA<coupledFvPatch>(fvp) ||
            isA<symmetryFvPatch>(fvp) ||
            isA<symmetryPlaneFvPatch>(fvp)
        )
        {
            notConstrain = false;
        }

        if (notConstrain)
        {
            gradIVF.boundaryFieldRef()[ipatch] = nf_.boundaryField()[ipatch]*iVF.boundaryField()[ipatch].snGrad();
        }
    }

    if(!Pstream::parRun())
    {
        return tgradIVF;
    }

    List<List3<scalar>> procVfValues(nProcPatches_);

    formVfValues(iVF,procVfValues);

    faceScalarDer(procVfValues[0],sVF.component(0),0,tField);
    gradIVF.boundaryFieldRef().replace(0,tField.boundaryField());

    faceScalarDer(procVfValues[0],sVF.component(0),1,tField);
    gradIVF.boundaryFieldRef().replace(3,tField.boundaryField());

    faceScalarDer(procVfValues[0],sVF.component(0),2,tField);
    gradIVF.boundaryFieldRef().replace(6,tField.boundaryField());

    faceScalarDer(procVfValues[1],sVF.component(1),0,tField);
    gradIVF.boundaryFieldRef().replace(1,tField.boundaryField());

    faceScalarDer(procVfValues[1],sVF.component(1),1,tField);
    gradIVF.boundaryFieldRef().replace(4,tField.boundaryField());

    faceScalarDer(procVfValues[1],sVF.component(1),2,tField);
    gradIVF.boundaryFieldRef().replace(7,tField.boundaryField());

    faceScalarDer(procVfValues[2],sVF.component(2),0,tField);
    gradIVF.boundaryFieldRef().replace(2,tField.boundaryField());

    faceScalarDer(procVfValues[2],sVF.component(2),1,tField);
    gradIVF.boundaryFieldRef().replace(5,tField.boundaryField());

    faceScalarDer(procVfValues[2],sVF.component(2),2,tField);
    gradIVF.boundaryFieldRef().replace(8,tField.boundaryField());

    return tgradIVF;
    
};

//- Calculate divergence of volume vector field on the faces.
//
// \param iVF        Internal vector field.
//                   Allowable values: constant reference to the volVectorField.
//
// \return           Divergence of iVF (scalar field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceScalarField> Foam::fvsc::leastSquaresOpt::Div(const volVectorField& iVF)
{

    surfaceVectorField sVF = linearInterpolate(iVF);
    surfaceScalarField tField = sVF.component(0)*0;
    tmp<surfaceScalarField> tdivIVF(0*nf_&fvc::snGrad(iVF));
    surfaceScalarField& divIVF = tdivIVF.ref();

    faceScalarDer(iVF.primitiveField().component(0),sVF.primitiveField().component(0),0,tField);
    divIVF.primitiveFieldRef() = tField;

    faceScalarDer(iVF.primitiveField().component(1),sVF.primitiveField().component(1),1,tField);
    divIVF.primitiveFieldRef() += tField;

    faceScalarDer(iVF.primitiveField().component(2),sVF.primitiveField().component(2),2,tField);
    divIVF.primitiveFieldRef() += tField;

    forAll(mesh_.boundaryMesh(), ipatch)
    {
        bool notConstrain = true;
        const fvPatch& fvp = mesh_.boundary()[ipatch];
        if
        (
            isA<emptyFvPatch>(fvp) ||
            isA<wedgeFvPatch>(fvp) ||
            isA<coupledFvPatch>(fvp) ||
            isA<symmetryFvPatch>(fvp) ||
            isA<symmetryPlaneFvPatch>(fvp)
        )
        {
            notConstrain = false;
        }

        if (notConstrain)
        {
            divIVF.boundaryFieldRef()[ipatch] = nf_.boundaryField()[ipatch] & iVF.boundaryField()[ipatch].snGrad();
        }
    }

    if(!Pstream::parRun())
    {
        return tdivIVF;
    }

    List<List3<scalar>> procVfValues(nProcPatches_); 
    formVfValues(iVF,procVfValues);

    
    faceScalarDer(procVfValues[0],sVF.component(0),0,tField);
    divIVF.boundaryFieldRef() = tField.boundaryField();

    faceScalarDer(procVfValues[1],sVF.component(1),1,tField);
    divIVF.boundaryFieldRef() += tField.boundaryField();

    faceScalarDer(procVfValues[2],sVF.component(2),2,tField);
    divIVF.boundaryFieldRef() += tField.boundaryField();

    
    return tdivIVF;
};

//- Calculate divergence of volume tensor field on the faces.
//
// \param iTF        Internal tensor field.
//                   Allowable values: constant reference to the volTensorField.
//
// \return           Divergence of iTF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::leastSquaresOpt::Div(const volTensorField& iTF)
{

    surfaceTensorField sTF = linearInterpolate(iTF);
    surfaceScalarField tField = sTF.component(0)*0;
    tmp<surfaceVectorField> tdivITF(0*nf_*fvc::snGrad(iTF.component(0)));
    surfaceVectorField& divITF = tdivITF.ref();
    surfaceScalarField divComp = tField;

    faceScalarDer(iTF.primitiveField().component(0),sTF.primitiveField().component(0),0,tField);
    divComp = tField;
    faceScalarDer(iTF.primitiveField().component(1),sTF.primitiveField().component(1),1,tField);
    divComp += tField;
    faceScalarDer(iTF.primitiveField().component(2),sTF.primitiveField().component(2),2,tField);
    divComp += tField;
    divITF.primitiveFieldRef().replace(0,divComp);

    faceScalarDer(iTF.primitiveField().component(3),sTF.primitiveField().component(3),0,tField);
    divComp = tField;
    faceScalarDer(iTF.primitiveField().component(4),sTF.primitiveField().component(4),1,tField);
    divComp += tField;
    faceScalarDer(iTF.primitiveField().component(5),sTF.primitiveField().component(5),2,tField);
    divComp += tField;
    divITF.primitiveFieldRef().replace(1,divComp);

    faceScalarDer(iTF.primitiveField().component(6),sTF.primitiveField().component(6),0,tField);
    divComp = tField;
    faceScalarDer(iTF.primitiveField().component(7),sTF.primitiveField().component(7),1,tField);
    divComp += tField;
    faceScalarDer(iTF.primitiveField().component(8),sTF.primitiveField().component(8),2,tField);
    divComp += tField;
    divITF.primitiveFieldRef().replace(2,divComp);

    forAll(mesh_.boundaryMesh(), ipatch)
    {
        bool notConstrain = true;
        const fvPatch& fvp = mesh_.boundary()[ipatch];
        if
        (
            isA<emptyFvPatch>(fvp) ||
            isA<wedgeFvPatch>(fvp) ||
            isA<coupledFvPatch>(fvp) ||
            isA<symmetryFvPatch>(fvp) ||
            isA<symmetryPlaneFvPatch>(fvp)
        )
        {
            notConstrain = false;
        }

        if (notConstrain)
        {
            divITF.boundaryFieldRef()[ipatch] = nf_.boundaryField()[ipatch] 
            & iTF.boundaryField()[ipatch].snGrad();
        }
    }


    if (!Pstream::parRun()) 
    {
        return tdivITF;
    }


    List<List3<scalar>> procVfValues(nProcPatches_);
    formVfValues(iTF,procVfValues);

    faceScalarDer(procVfValues[0],sTF.component(0),0,tField);
    divComp = tField;
    faceScalarDer(procVfValues[1],sTF.component(1),1,tField);
    divComp += tField;
    faceScalarDer(procVfValues[2],sTF.component(2),2,tField);
    divComp += tField;
    divITF.boundaryFieldRef().replace(0,divComp.boundaryField());

    faceScalarDer(procVfValues[3],sTF.component(3),0,tField);
    divComp = tField;
    faceScalarDer(procVfValues[4],sTF.component(4),1,tField);
    divComp += tField;
    faceScalarDer(procVfValues[5],sTF.component(5),2,tField);
    divComp += tField;
    divITF.boundaryFieldRef().replace(1,divComp.boundaryField());

    faceScalarDer(procVfValues[6],sTF.component(6),0,tField);
    divComp = tField;
    faceScalarDer(procVfValues[7],sTF.component(7),1,tField);
    divComp += tField;
    faceScalarDer(procVfValues[8],sTF.component(8),2,tField);
    divComp += tField;
    divITF.boundaryFieldRef().replace(2,divComp.boundaryField());

    return tdivITF;
}

//
//END-OF-FILE
//


