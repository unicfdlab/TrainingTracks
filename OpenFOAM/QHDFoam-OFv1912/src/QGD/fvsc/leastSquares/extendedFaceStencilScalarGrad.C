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
    Foam::fvsc::leastSquares::extendedFaceStencilScalarGrad
\*---------------------------------------------------------------------------*/

#include "leastSquaresStencil.H"
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
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::leastSquares::Grad(const volScalarField& iF)
{
    surfaceScalarField sF = linearInterpolate(iF);
    surfaceScalarField sngF (fvc::snGrad(iF));

    tmp<surfaceVectorField> tgradIF(0.0 * nf_ * sngF);
    surfaceVectorField& gradIF = tgradIF.ref();
    
    // List of faces
    const faceList& faces = mesh_.faces();

    vector gf = vector::zero;
    forAll(faces, facei)
    {
        if (mesh_.isInternalFace(facei))
        {
            gf = vector::zero;
            forAll(GdfAll_[facei], i)
            {
                gf = gf + wf2All_[facei][i]*GdfAll_[facei][i]*(iF[neighbourCells_[facei][i]] - sF[facei]);
            }

            gradIF[facei] = gf;
        }
    }
    
    //process faces with degenerate stencil
    label dFaceId = -1;
    forAll(internalDegFaces_, facei)
    {
        dFaceId = internalDegFaces_[facei];
        //gradIF[dFaceId] = nf_[dFaceId] * sngF[dFaceId];
        gradIF[dFaceId] = sngF[dFaceId] * nf_[dFaceId];
    }
    
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
//            fvp.coupled()
        )
        {
            notConstrain = false;
        }

        if (notConstrain)
        {
            gradIF.boundaryFieldRef()[ipatch] = nf_.boundaryField()[ipatch] * 
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
    //set values from this domain
    label cellId = -1;
    forAll(procPairs_, patchI)
    {
        if (procPairs_[patchI] > -1)
        {
            procVfValues[patchI].resize(procWf2_[patchI].size());
            forAll(procVfValues[patchI], faceI)
            {
                procVfValues[patchI][faceI].resize(procWf2_[patchI][faceI].size());
                procVfValues[patchI][faceI] = 0.0; //make values zero
                forAll(myProcPatchCells_[patchI][faceI], cellI)
                {
                    cellId = myProcPatchCells_[patchI][faceI][cellI];
                    procVfValues[patchI][faceI][cellI] = iF.primitiveField()[cellId];
                }
            }
        }
    }
    
    //Step 1. Send field data to neighbouring processors (non-blocking mode)
    
    PstreamBuffers pBuffers(Pstream::commsTypes::nonBlocking);
    forAll(procPairs_, procI)
    {
        label procId = neigProcs_[procI];
//        label dataSz = 0;
        
        DynamicList<scalar> locVf;
        
        if (procPairs_[procI] > -1) //patch proc pair
        {
            forAll(procVfValues[procI], faceI)
            {
                for(
                        label
                        cellI = 0;
                        cellI <= ownEnd_[procI][faceI];
                        cellI++
                    )
                {
                    locVf.append(procVfValues[procI][faceI][cellI]);
                }

            }
        }
        else //corner connected process
        {
            label cellId = -1;
            label addrId = corProcIds_[procId];
            //label addrId = corProcIds_.capacity();
            forAll(corCellIds_[addrId], iCellId)
            {
                cellId = corCellIds_[addrId][iCellId];
                locVf.append(iF.primitiveField()[cellId]);
            }
        }
        
        UOPstream oProcStr(procId, pBuffers);
        oProcStr << locVf;
    }
    
    //Step 2. Recieve field data from neighbouring processors
    pBuffers.finishedSends();
    label iCorProc = 0;
    forAll(procPairs_, procI)
    {
        label procId = neigProcs_[procI];
        
        UIPstream iProcStr(procId, pBuffers);
        List<scalar> locVf (iProcStr);
        
        if (procPairs_[procI] > -1)
        {
            label iVf = 0;
            forAll(neiStart_[procI], iFace)
            {
                for(
                        label
                        iCell=neiStart_[procI][iFace];
                        iCell<=neiEnd_[procI][iFace];
                        iCell++
                    )
                {
                    procVfValues[procI][iFace][iCell] = 
                        locVf[iVf];
                    iVf++;
                }
            }
        }
        else
        {
            label patchNo = -1;
            label faceNo  = -1;
            label cellNo  = -1;
            label offset  = -1;
            
            const List<Triple<label> >& addr = corAddr_[iCorProc];
            
            forAll(addr, iVal)
            {
                patchNo = addr[iVal][0];
                faceNo  = addr[iVal][1];
                cellNo  = addr[iVal][2];
                
                offset = corStart_[patchNo][faceNo];
                procVfValues[patchNo][faceNo][cellNo+offset] = locVf[iVal];
            }
            iCorProc++;
        }
    }

    //Step 3. Calculate gradient at faces on processor patches
    forAll(procPairs_, patchI)
    {
        label procPatchId = procPairs_[patchI];
        if (procPatchId > -1)
        {
            fvsPatchVectorField& pgradf = gradIF.boundaryFieldRef()[procPatchId];
            const List2<scalar> & pvf = procVfValues[patchI];
            const List2<scalar> & pwf2= procWf2_[patchI];
            const List2<vector> & pgdf= procGdf_[patchI];
        
            forAll(pgradf, iFace)
            {
                gf = vector::zero;
                forAll(procGdf_[patchI][iFace], i)
                {
                    gf = gf + pwf2[iFace][i]*pgdf[iFace][i]*
                            (pvf[iFace][i] - sF.boundaryField()[procPatchId][iFace]);
                }
                pgradf[iFace] = gf;
            }
            
            //Update processor degenerate faces
            const labelList& degProcFaces = procDegFaces_[patchI];
            label degId = -1;
            forAll(degProcFaces, iFace)
            {
                degId = degProcFaces[iFace];
                
                pgradf[degId] = nf_.boundaryField()[procPatchId][degId]*
                    sngF.boundaryField()[procPatchId][degId];
            }
        }
    }
    
    return tgradIF;
};

//
//END-OF-FILE
//


