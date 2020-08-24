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
    grpleastSquares
Class
   Foam::fvsc::leastSquares::leastSquaresStencil
Description 
   Realisation least squares method for calculationg of differential operators
\*---------------------------------------------------------------------------*/

#include "leastSquaresStencil.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <HashTable.H>
#include "addToRunTimeSelectionTable.H"
#include "faceSet.H"

namespace Foam
{
namespace fvsc
{
    defineTypeNameAndDebug(leastSquares,0);
    addToRunTimeSelectionTable
    (
        fvscStencil,
        leastSquares,
        components
    );
}
}


// constructors
Foam::fvsc::leastSquares::leastSquares(const IOobject& io)
:
    fvscStencil(io),
    leastSquaresBase(mesh_)
{
    faceSet degenerateFacesSet
    (
        //faceSetHeader
        mesh_,
        "degenerateStencilFaces",
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );
    
    if (degenerateFacesSet.size() > 0) //.headerOk())
    {
        Info << "Found set with faces for reduced approximation QGD terms" << endl;
        
        //read list of degenerated faces

    
        const labelList degenerateFaces = degenerateFacesSet.toc();
        
        labelHashSet intDegFaces;
        List<labelHashSet> procDegFaces (procDegFaces_.size());
    
        forAll(degenerateFaces, iDegFace)
        {
            label faceId = degenerateFaces[iDegFace];
    
            if (mesh_.isInternalFace(faceId))
            {
                intDegFaces.insert(faceId);
            }
            else
            {
                label patchId = mesh_.boundaryMesh().whichPatch(faceId);
                label iPatch = -1;
                forAll(procPairs_, iPP)
                {
                    if (procPairs_[iPP] == patchId)
                    {
                        iPatch = iPP;
                        break;
                    }
                }
                
                if (iPatch > -1)
                {
                    procDegFaces[iPatch].insert
                    (
                        faceId
                        -
                        mesh_.boundaryMesh()[patchId].start()
                    );
                }
            }
        }
        
        internalDegFaces_.append(intDegFaces.toc());
        forAll(procDegFaces, iProcPair)
        {
            if (procDegFaces[iProcPair].toc().size() > 0)
            {
                procDegFaces_[iProcPair].append
                (
                    procDegFaces[iProcPair].toc()
                );
            }
        }
    }
    else
    {
        Info << "Set \"degenerateStencilFaces\" with faces for reduced approximation QGD terms not found" << endl;
    }
}

Foam::fvsc::leastSquares::~leastSquares()
{
}

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceTensorField> Foam::fvsc::leastSquares::Grad(const volVectorField& iVF)
{
    surfaceVectorField gradComp0col = Grad(iVF.component(0));
    surfaceVectorField gradComp1col = Grad(iVF.component(1));
    surfaceVectorField gradComp2col = Grad(iVF.component(2));

    tmp<surfaceTensorField> tgradIVF(0* nf_* fvc::snGrad(iVF));
    surfaceTensorField& gradIVF = tgradIVF.ref();
    
    //set internal field
    gradIVF.primitiveFieldRef().replace(0, gradComp0col.primitiveField().component(0));
    gradIVF.primitiveFieldRef().replace(1, gradComp1col.primitiveField().component(0));
    gradIVF.primitiveFieldRef().replace(2, gradComp2col.primitiveField().component(0));
    
    gradIVF.primitiveFieldRef().replace(3, gradComp0col.primitiveField().component(1));
    gradIVF.primitiveFieldRef().replace(4, gradComp1col.primitiveField().component(1));
    gradIVF.primitiveFieldRef().replace(5, gradComp2col.primitiveField().component(1));
    
    gradIVF.primitiveFieldRef().replace(6, gradComp0col.primitiveField().component(2));
    gradIVF.primitiveFieldRef().replace(7, gradComp1col.primitiveField().component(2));
    gradIVF.primitiveFieldRef().replace(8, gradComp2col.primitiveField().component(2));
    
    //set external fields
    forAll(mesh_.boundaryMesh(), patchi)
    {
        forAll(mesh_.boundary()[patchi], facei)
        {
            gradIVF.boundaryFieldRef()[patchi][facei].component(0) = 
                gradComp0col.boundaryField()[patchi][facei].component(0);
            gradIVF.boundaryFieldRef()[patchi][facei].component(1) = 
                gradComp1col.boundaryField()[patchi][facei].component(0);
            gradIVF.boundaryFieldRef()[patchi][facei].component(2) = 
                gradComp2col.boundaryField()[patchi][facei].component(0);

            gradIVF.boundaryFieldRef()[patchi][facei].component(3) = 
                gradComp0col.boundaryField()[patchi][facei].component(1);
            gradIVF.boundaryFieldRef()[patchi][facei].component(4) = 
                gradComp1col.boundaryField()[patchi][facei].component(1);
            gradIVF.boundaryFieldRef()[patchi][facei].component(5) = 
                    gradComp2col.boundaryField()[patchi][facei].component(1);

            gradIVF.boundaryFieldRef()[patchi][facei].component(6) = 
                gradComp0col.boundaryField()[patchi][facei].component(2);
            gradIVF.boundaryFieldRef()[patchi][facei].component(7) = 
                gradComp1col.boundaryField()[patchi][facei].component(2);
            gradIVF.boundaryFieldRef()[patchi][facei].component(8) = 
                gradComp2col.boundaryField()[patchi][facei].component(2);
        }
    }

    return tgradIVF;
};

//- Calculate divergence of volume vector field on the faces.
//
// \param iVF        Internal vector field.
//                   Allowable values: constant reference to the volVectorField.
//
// \return           Divergence of iVF (scalar field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceScalarField> Foam::fvsc::leastSquares::Div(const volVectorField& iVF)
{
    surfaceVectorField gradComp0 = Grad(iVF.component(0));
    surfaceVectorField gradComp1 = Grad(iVF.component(1));
    surfaceVectorField gradComp2 = Grad(iVF.component(2));

    tmp<surfaceScalarField> tdivIVF(0 * (nf_ & fvc::snGrad(iVF)));
    surfaceScalarField& divIVF = tdivIVF.ref();
    
    divIVF.primitiveFieldRef() = gradComp0.primitiveField().component(0)
                               + gradComp1.primitiveField().component(1)
                               + gradComp2.primitiveField().component(2);
    
    forAll(mesh_.boundary(), patchi)
    {
        divIVF.boundaryFieldRef()[patchi] = 
            gradComp0.boundaryField()[patchi].component(0)
            +
            gradComp1.boundaryField()[patchi].component(1)
            +
            gradComp2.boundaryField()[patchi].component(2);
    }
    
    return tdivIVF;
};

//- Calculate divergence of volume tensor field on the faces.
//
// \param iTF        Internal tensor field.
//                   Allowable values: constant reference to the volTensorField.
//
// \return           Divergence of iTF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::leastSquares::Div(const volTensorField& iTF)
{
    tmp<surfaceVectorField> gradComp0 (Grad(iTF.component(0)));
    tmp<surfaceVectorField> gradComp1 (Grad(iTF.component(1)));
    tmp<surfaceVectorField> gradComp2 (Grad(iTF.component(2)));
    
    tmp<surfaceVectorField> gradComp3 (Grad(iTF.component(3)));
    tmp<surfaceVectorField> gradComp4 (Grad(iTF.component(4)));
    tmp<surfaceVectorField> gradComp5 (Grad(iTF.component(5)));

    tmp<surfaceVectorField> gradComp6 (Grad(iTF.component(6)));
    tmp<surfaceVectorField> gradComp7 (Grad(iTF.component(7)));
    tmp<surfaceVectorField> gradComp8 (Grad(iTF.component(8)));

    tmp<surfaceScalarField> divComp0 (gradComp0().component(0) + gradComp3().component(1) + gradComp6().component(2));
    tmp<surfaceScalarField> divComp1 (gradComp1().component(0) + gradComp4().component(1) + gradComp7().component(2));
    tmp<surfaceScalarField> divComp2 (gradComp2().component(0) + gradComp5().component(1) + gradComp8().component(2));

    tmp<surfaceVectorField> tdivITF(0*nf_*fvc::snGrad(iTF.component(0)));
    surfaceVectorField& divITF = tdivITF.ref();
    
    divITF.primitiveFieldRef().replace(0, divComp0().primitiveField());
    divITF.primitiveFieldRef().replace(1, divComp1().primitiveField());
    divITF.primitiveFieldRef().replace(2, divComp2().primitiveField());
    
    forAll(mesh_.boundary(), patchi)
    {
        forAll(mesh_.boundary()[patchi], facei)
        {
            divITF.boundaryFieldRef()[patchi][facei].component(0) = 
                divComp0().boundaryField()[patchi][facei];
            divITF.boundaryFieldRef()[patchi][facei].component(1) = 
                divComp1().boundaryField()[patchi][facei];
            divITF.boundaryFieldRef()[patchi][facei].component(2) = 
                divComp2().boundaryField()[patchi][facei];
        }
    }
    
    return tdivITF;
}

//
//END-OF-FILE
//


