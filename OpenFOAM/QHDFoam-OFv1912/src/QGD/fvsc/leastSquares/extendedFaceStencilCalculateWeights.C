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
    Foam::fvsc::leastSquares::extendedFaceStencilCalculateWeights
Description 
    Base methods for calculating weights
\*---------------------------------------------------------------------------*/



#include "leastSquaresBase.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <HashTable.H>

//- Compute weights for least squares scheme for gradient calculation.
void Foam::fvsc::leastSquaresBase::calculateWeights()
{
    //Pout << "Start calculateWeights()" << endl;
    
    const faceList& faces = cMesh_.faces();
    GdfAll_.resize(faces.size());
    wf2All_.resize(faces.size());
    label cellDim = 3;
    
    //create list of underdetermined cells:
    //cellSet udCells(mesh_, "udCells", mesh_.nCells()/100);
    //mesh_.checkCellDeterminant(true, &udCells);
    label nDegFaces = 0;
    scalar minDet = GREAT;
    scalar detG = 0.0;
    
    //for internal faces
    forAll(faces, facei)
    {
        if (cMesh_.isInternalFace(facei))
        {
            List<vector> df(neighbourCells_[facei].size());
            scalarList wf2(neighbourCells_[facei].size());
            symmTensor G(0);
            
            vector Cf = cMesh_.faceCentres()[facei];
            
            forAll(neighbourCells_[facei], i)
            {
                df[i] = cMesh_.cellCentres()[neighbourCells_[facei][i]] - Cf;
                wf2[i] = 1/magSqr(df[i]);
                symmTensor addToG(0);
                addToG = sqr(df[i]);
                addToG = addToG * wf2[i];
                G += addToG;
            }
            
            symmTensor G0(0);
            cellDim = 3;
            
            //correct G tensor for zero directions
            {
                if (mag(G.xx()) < SMALL)
                {
                    G0.xx() = 1;
                    cellDim--;
                }
                if (mag(G.yy()) < SMALL)
                {
                    G0.yy() = 1;
                    cellDim--;
                }
                if (mag(G.zz()) < SMALL)
                {
                    G0.zz() = 1;
                    cellDim--;
                }
                
                if (cellDim != cMesh_.nGeometricD())
                {
                    WarningInFunction
                        << "face " << facei << " with center "
                        << cMesh_.faceCentres()[facei] << nl
                        << " connected to cells with dimensions " << cellDim
                        << " less then geometric " << cMesh_.nGeometricD() 
                        << nl << endl;
                    Pout << "Degenerate face: " << facei << endl;
                }
            }
            
            /*
            if(mesh_.nGeometricD()==1)
            {
                symmTensor G01(1, 0, 0, 1, 0, 1);
                symmTensor G02(sqr((Vector<label>::one + mesh_.geometricD())/2));
                G0 = G01 - G02;
            }
            else
            {
                if(mesh_.nGeometricD()==2)
                {
                    G0 = sqr((Vector<label>::one - mesh_.geometricD())/2);
                }
            };
            */

            G = G + G0;
            detG = det(G);
            if (detG < minDet)
            {
                minDet = detG;
            }
            
            if (detG < 1)
            {
                nDegFaces++;
                internalDegFaces_.append(facei);
            }
            else
            {
                G = inv(G);
                G = G - G0;
            }
            
            forAll(df, i)
            {
                df[i] = G&df[i];
            }

            GdfAll_[facei] = df;
            wf2All_[facei] = wf2;
        }
    }
    
    if (nDegFaces > 0)
    {
        Pout << "Min determinant      : " << minDet << endl;
        Pout << "Total # of deg. faces: " << nDegFaces << endl;
    }
    
    //Info << "End for not parallel" << endl;
    
    if (Pstream::parRun())
    {
        List3<vector> ownCellCenters(nProcPatches_);
        List3<vector> neiCellCenters(nProcPatches_);
        List3<vector> corCellCenters(nProcPatches_);
        //List<label> nOwnCells(nProcPatches_, 0);
        
        label cellId = -1;
        label nCorCells = -1;
        
        //Step 1. Set cell centers at my patches
        //Cell centers are stored only per processor patch
        forAll(ownCellCenters, iProcPatch)
        {
            ownCellCenters[iProcPatch].resize(myProcPatchCells_[iProcPatch].size());
            neiCellCenters[iProcPatch].resize(myProcPatchCells_[iProcPatch].size());
            corCellCenters[iProcPatch].resize(myProcPatchCells_[iProcPatch].size());
            forAll(ownCellCenters[iProcPatch], iFace)
            {
                ownCellCenters[iProcPatch][iFace].resize
                (
                    myProcPatchCells_[iProcPatch][iFace].size()
                );
                
                nCorCells = corEnd_[iProcPatch][iFace] - corStart_[iProcPatch][iFace] + 1;
                
                //Pout << "corEnd = " << corEnd_[iProcPatch][iFace]  << " nCorCell = " << nCorCells << endl;
                corCellCenters[iProcPatch][iFace].resize
                (
                    nCorCells
                );
                
                forAll(ownCellCenters[iProcPatch][iFace], iCell)
                {
                    cellId = myProcPatchCells_[iProcPatch][iFace][iCell];
                    ownCellCenters[iProcPatch][iFace][iCell] = cMesh_.C()[cellId];
                }
                    //nOwnCells[iProcPatch] += ownCellCenters[iProcPatch][iFace].size();
            }
        }

        // Step 2. Loop over all neighboring processors and send/receive cell centers
        {
            PstreamBuffers pBuffers(Pstream::commsTypes::nonBlocking);
        
            forAll(neigProcs_, iProcPair)
            {
                if (procPairs_[iProcPair] > -1) //send cell centers for patch-neighbouring processes
                {
                    label procId = neigProcs_[iProcPair];
                    UOPstream oProcStr(procId, pBuffers);
                    oProcStr << ownCellCenters[iProcPair];
                }
            }
            
            pBuffers.finishedSends();
        
            forAll(neigProcs_, iProcPair)
            {
                if (procPairs_[iProcPair] > -1) //recieve cell centers for patch-neighbouring processes
                {
                    label procId = neigProcs_[iProcPair];
                    UIPstream iProcStr(procId, pBuffers);
                    iProcStr >> neiCellCenters[iProcPair];
                }
            }
            
            //Pout << "neiCellCenters = " << neiCellCenters << endl;
        }
        
        // Step 3. Loop over all corner neigbouring processors and send/receive cell centers
        {
            PstreamBuffers pBuffers(Pstream::commsTypes::nonBlocking);
            
            // Send
            forAll(neigProcs_, iProcPair)
            {
                if (procPairs_[iProcPair] < 0)
                {
                    label procId = neigProcs_[iProcPair];
                    UOPstream oProcStr(procId, pBuffers);
                    label id = corProcIds_[procId];
                    
                    List<vector> locCc (corCellIds_[id].size());
                    
                    forAll(locCc, iCell)
                    {
                        label cellId = corCellIds_[id][iCell];
                        locCc[iCell] = cMesh_.C()[cellId];
                    }
                    oProcStr << locCc;
                    //Pout << "Sending " << locCc << " to " << procId << endl;
                }
            }
            
            // Recieve
            pBuffers.finishedSends();
            label iCorProc = 0;
            forAll(neigProcs_, iProcPair)
            {
                if (procPairs_[iProcPair] < 0)
                {
                    label procId = neigProcs_[iProcPair];
                    UIPstream iProcStr(procId, pBuffers);
                    
                    List<vector> corCc (iProcStr);
                    
                    //Pout << "Received from " << procId << " cell centers " << corCc << endl;
                    
                    const List<Triple<label> > & addr = corAddr_[iCorProc];
                    label patchNo = -1;
                    label faceNo  = -1;
                    label cellNo  = -1;
                    
                    forAll(corCc, iCell)
                    {
                        patchNo = addr[iCell][0];
                        faceNo  = addr[iCell][1];
                        cellNo  = addr[iCell][2];
                        
                        corCellCenters[patchNo][faceNo][cellNo] = corCc[iCell];
                    }
                    
                    iCorProc++;
                }
            }
            
            //Pout << "corCellCenters = " << corCellCenters << endl;
        }
        
        
        // Step 4. Calculate weights
        //Pout << "neiCellCenters = " << neiCellCenters << endl;
        //Pout << "corCellCenters = " << corCellCenters << endl;
        procDegFaces_.resize(nProcPatches_);
        forAll(myProcPatchCells_,iProcPatch)
        {
            if (procPairs_[iProcPatch] > -1)
            {
                const label iProcPatchId = procPairs_[iProcPatch];
                const fvPatch& fvp       = cMesh_.boundary()[iProcPatchId];
                
                label nFaceCells = 0;
                
                forAll(ownCellCenters[iProcPatch], facei)
                {
                    nFaceCells = corEnd_[iProcPatch][facei] + 1;
                    
                    List<vector> df (nFaceCells, vector::zero);
                    List<scalar> wf2(nFaceCells, 0.0);

                    symmTensor G(0);
                    
                    forAll(ownCellCenters[iProcPatch][facei], i)
                    {
                        df[i] = ownCellCenters[iProcPatch][facei][i] - fvp.Cf()[facei];
                        wf2[i] = 1/magSqr(df[i]);
                        symmTensor addToG(0);
                        addToG = sqr(df[i]);
                        addToG = addToG * wf2[i];
                        G += addToG;
                    }
                    
                    label k = neiStart_[iProcPatch][facei];
                    forAll(neiCellCenters[iProcPatch][facei], i)
                    {
                        df[k] = neiCellCenters[iProcPatch][facei][i] - fvp.Cf()[facei];
                        wf2[k] = 1/magSqr(df[k]);
                        symmTensor addToG(0);
                        addToG = sqr(df[k]);
                        addToG = addToG * wf2[k];
                        G += addToG;
                        k++;
                    }
                    
                    label l = corStart_[iProcPatch][facei];
                    forAll(corCellCenters[iProcPatch][facei], i)
                    {
                        df[l] = corCellCenters[iProcPatch][facei][i] - fvp.Cf()[facei];
                        wf2[l] = 1/magSqr(df[l]);
                        symmTensor addToG(0);
                        addToG = sqr(df[l]);
                        addToG = addToG * wf2[l];
                        G += addToG;
                        l++;
                    }
                    
                    symmTensor G0(0);
                    cellDim = 3;
                    /*
                    if(mesh_.nGeometricD()==1)
                    {
                        symmTensor G01(1, 0, 0, 1, 0, 1);
                        symmTensor G02(sqr((Vector<label>::one + mesh_.geometricD())/2));
                        G0 = G01 - G02;
                    }
                    else
                    {
                        if(mesh_.nGeometricD()==2)
                        {
                            G0 = sqr((Vector<label>::one - mesh_.geometricD())/2);
                        }
                    };
                    */
        
                    //correct G tensor for zero directions
                    {
                        if (mag(G.xx()) < SMALL)
                        {
                            G0.xx() = 1;
                            cellDim--;
                        }
                        if (mag(G.yy()) < SMALL)
                        {
                            G0.yy() = 1;
                            cellDim--;
                        }
                        if (mag(G.zz()) < SMALL)
                        {
                            G0.zz() = 1;
                            cellDim--;
                        }
                    }

                    if (cellDim != cMesh_.nGeometricD())
                    {
                        WarningInFunction
                            << "face " << facei << " with center "
                            << fvp.Cf()[facei] << nl
                            << " connected to cells with dimensions " << cellDim
                            << " less then geometric " << cMesh_.nGeometricD() 
                            << nl << endl;
                    }
                    
                    G = G + G0;
                    
                    detG = det(G);
                    if (detG < 1)
                    {
                        procDegFaces_[iProcPatch].append(facei);
                    }
                    else
                    {
                        G = inv(G);
                        G = G - G0;
                    }
                    
                    forAll(df, i)
                    {
                        df[i] = G&df[i];
                    }
                    
                    procGdf_[iProcPatch][facei] = df;
                    procWf2_[iProcPatch][facei] = wf2;
                }
            }
        }
    }
    
    //Pout << "End calculateWeights()" << endl;
};

