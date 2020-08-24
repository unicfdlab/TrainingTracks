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
    Foam::fvsc::leastSquaresOpt::extendedFaceStencilScalarGradOpt
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
#include <HashTable.H>

#include "emptyFvPatch.H"
#include "coupledFvPatch.H"
#include "wedgeFvPatch.H"
#include "symmetryFvPatch.H"
#include "symmetryPlaneFvPatch.H"

//- Calculate gradient of volume scalar function on the faces
//
// \param iF         Internal scalar field.
//                   Allowable values: constant reference to the volScalarField.
//
// \return           Gradient of iF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::leastSquaresOpt::Grad(const volScalarField& iF)
{
    surfaceScalarField sF = linearInterpolate(iF); 
    return Grad(iF,sF);
};

Foam::tmp<Foam::surfaceVectorField> 
Foam::fvsc::leastSquaresOpt::Grad(const volScalarField& iF, const surfaceScalarField& sF)
{

    tmp<surfaceVectorField> tgradIF(0.0*nf_*fvc::snGrad(iF));
    surfaceVectorField& gradIF = tgradIF.ref();
    //scalarField tField = sF;
    surfaceScalarField tField = sF*0;

    faceScalarDer(iF.primitiveField(),sF.primitiveField(),0,tField);
    gradIF.primitiveFieldRef().replace(0, tField);
    faceScalarDer(iF.primitiveField(),sF.primitiveField(),1,tField);
    gradIF.primitiveFieldRef().replace(1, tField);
    faceScalarDer(iF.primitiveField(),sF.primitiveField(),2,tField);
    gradIF.primitiveFieldRef().replace(2, tField);
    
    //update boundary field
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
            gradIF.boundaryFieldRef()[ipatch] = 
                nf_.boundaryField()[ipatch] * 
                iF.boundaryField()[ipatch].snGrad();
        }
    }

    if(!Pstream::parRun())
    {
        return tgradIF;
    }
    
    /*
     *
     * Update processor patches for parallel case
     *
     */
    
    //allocate storage for near-patch field

    List3<scalar> procVfValues(nProcPatches_); //array of values from neighb. processors
    formVfValues(iF,procVfValues);

    faceScalarDer(procVfValues,sF,0,tField);
    gradIF.boundaryFieldRef().replace(0, tField.boundaryFieldRef());
    faceScalarDer(procVfValues,sF,1,tField);
    gradIF.boundaryFieldRef().replace(1, tField.boundaryFieldRef());
    faceScalarDer(procVfValues,sF,2,tField);
    gradIF.boundaryFieldRef().replace(2, tField.boundaryFieldRef());
    
    return tgradIF;


};

//
//END-OF-FILE
//

