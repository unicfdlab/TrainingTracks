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
    Foam::fvsc::leastSquares::extendedFaceStencilFindNeighbours
Description 
    Base methods for finding neighbours
\*---------------------------------------------------------------------------*/
#include "leastSquaresBase.H"
#include "processorFvPatch.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <HashTable.H>

//- Find neighbour cells for each face (throught face points).
void Foam::fvsc::leastSquaresBase::findNeighbours()
{
    // List of faces
    const faceList& faces = cMesh_.faces();
    labelListList neighbourCellsForFace(cMesh_.nInternalFaces());
    
    //Step 1. Find neighbours for internal faces
    forAll(faces, facei)
    {
        if(cMesh_.isInternalFace(facei))
        {
            labelList neighbourCells;
            labelList pointsFacei = faces[facei];

            forAll(pointsFacei, k)
            {
                label pointi = pointsFacei[k];
                // cells should be added to list if it wasn't contain in the list yet
                labelList neighbourCellsForPointI = cMesh_.pointCells()[pointi];

                forAll(neighbourCellsForPointI, j)
                {
                    label celli = neighbourCellsForPointI[j];

                    bool contained = false;

                    forAll(neighbourCells, f)
                    {
                        label ncell=neighbourCells[f];
                        if(ncell==celli)
                        {
                            contained = true;
                        }
                    }

                    if(!contained)
                    {
                        neighbourCells.append(celli);
                    }
                }
            }
            neighbourCellsForFace[facei].append(neighbourCells);
        }
    }
    
    neighbourCells_ = neighbourCellsForFace;
    
    //Step 2. Find neighbours for processors faces
    if (Pstream::parRun())
    {
        //Step 2.1 Load processor addresing
        labelIOList localPointProcAddr
        (
           IOobject
            (
                "pointProcAddressing",
                cMesh_.facesInstance(),
                cMesh_.meshSubDir,
                cMesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        labelHashTable<label> globalPointProcAddr; //reverse g==>l point addr
        forAll(localPointProcAddr, iPoint)
        {
            globalPointProcAddr.insert
            (
                localPointProcAddr[iPoint],
                iPoint
            );
        }
        
        labelIOList localCellProcAddr
        (
           IOobject
            (
                "cellProcAddressing",
                cMesh_.facesInstance(),
                cMesh_.meshSubDir,
                cMesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        labelHashTable<label> globalCellProcAddr; //reverse g==>l cell addr
        forAll(localCellProcAddr, iCell)
        {
            globalCellProcAddr.insert
            (
                localCellProcAddr[iCell],
                iCell
            );
        }
        
        forAll(cMesh_.boundary(), patchIndex)
        {
            if (isType<processorFvPatch>(cMesh_.boundary()[patchIndex]))
            {
                procPairs_.append(patchIndex);
                idProcPatchPairs_.insert
                (
                    patchIndex,
                    procPairs_.size() - 1
                );
            }
        }
        nProcPatches_ = procPairs_.size();
        procGdf_.resize(nProcPatches_);
        procWf2_.resize(nProcPatches_);
        neigProcs_.resize(nProcPatches_);
        
        labelHashTable<labelHashSet> locPointProcs;
        forAll(procPairs_, iProcPair)
        {
            const label iProcPatchId = procPairs_[iProcPair];
            const fvPatch& fvp       = cMesh_.boundary()[iProcPatchId];
            const processorFvPatch& procp = refCast<const processorFvPatch>(fvp);
            procGdf_[iProcPair].resize(fvp.size());
            procWf2_[iProcPair].resize(fvp.size());
            neigProcs_[iProcPair] = procp.neighbProcNo();
        }
        
        myProcPatchCells_.resize(nProcPatches_);
        forAll(myProcPatchCells_, iProcPatch)
        {
            //looking for cells connected to procesorPatch
            if (procPairs_[iProcPatch] > -1)
            {
                const label iProcPatchId = procPairs_[iProcPatch];
                const fvPatch& fvp       = cMesh_.boundary()[iProcPatchId];
                const polyPatch& pp      = fvp.patch(); //list of faces from which the patch is build up
                myProcPatchCells_[iProcPatch].resize(fvp.size());
                
                label pointId = -1;
                label gPointId = -1;
                forAll(pp, facei)
                {
                    labelList& faceCells = myProcPatchCells_[iProcPatch][facei];
                    labelHashSet vCells;
                    forAll(pp[facei], pointi)
                    {
                        pointId = pp[facei][pointi];
                        //add only unique cells ids
                        vCells.insert(cMesh_.pointCells()[pointId]);
                        //take global point id
                        gPointId = localPointProcAddr[pointId];
                        //
                        if (locPointProcs.found(gPointId))
                        {
                            locPointProcs[gPointId].insert(neigProcs_[iProcPatch]);
                        }
                        else
                        {
                            locPointProcs.insert(gPointId, labelHashSet());
                            locPointProcs[gPointId].insert(neigProcs_[iProcPatch]);
                        }
                    }
                    faceCells.append(vCells.toc());
                }
            }
        }
        
        //create and distribute internal addressing
        ownEnd_.resize(nProcPatches_);
        neiStart_.resize(nProcPatches_);
        neiEnd_.resize(nProcPatches_);
        corStart_.resize(nProcPatches_);
        corEnd_.resize(nProcPatches_);

        //Pout << "Addressing for patch neigbours created" << endl;
        
        {
            PstreamBuffers pBuffers (Pstream::commsTypes::nonBlocking);
            
            //Set and send
            forAll(myProcPatchCells_, iProcPatch)
            {
                if (procPairs_[iProcPatch] > -1)
                {
                    label nFaces = myProcPatchCells_[iProcPatch].size();
                    ownEnd_[iProcPatch].resize(nFaces);
                    neiStart_[iProcPatch].resize(nFaces);
                    forAll(myProcPatchCells_[iProcPatch], iFace)
                    {
                        ownEnd_[iProcPatch][iFace] = myProcPatchCells_[iProcPatch][iFace].size() - 1;
                        neiStart_[iProcPatch][iFace] = ownEnd_[iProcPatch][iFace] + 1;
                    }
                    
                    UOPstream oProcStr(neigProcs_[iProcPatch], pBuffers);
                    oProcStr << neiStart_[iProcPatch];
                }
            }
            
            pBuffers.finishedSends();
            
            //Recieve
            forAll(myProcPatchCells_, iProcPatch)
            {
                if (procPairs_[iProcPatch] > -1)
                {
                    label nFaces = myProcPatchCells_[iProcPatch].size();
                    List<label> neiLen (nFaces, 0);
                    neiEnd_[iProcPatch].resize(nFaces);
                    corStart_[iProcPatch].resize(nFaces);
                    corEnd_[iProcPatch].resize(nFaces);
                
                    UIPstream iProcStr(neigProcs_[iProcPatch], pBuffers);
                    iProcStr >> neiLen;
                
                    neiEnd_[iProcPatch] = ownEnd_[iProcPatch] + neiLen;
                    corStart_[iProcPatch] = neiEnd_[iProcPatch] + 1;
                    corEnd_[iProcPatch] = corStart_[iProcPatch] - 1;
                }
            }
        }

        //select point which belongs to more than 2 processors
        
        DynamicList<label> multipleProcsPoints;
        
        forAllConstIter(labelHashTable<labelHashSet>, locPointProcs, iter)
        {
            if (iter().size() > 1) //point at multiple processor's boundaries
            {
                multipleProcsPoints.append(iter.key());
            }
        }
        
        //accumulate cell id's for each corner point
        //at master process
        //fill cell-to-processor table
        {
            label gPointId = -1;
            label lPointId = -1;
            List<label> gPointCells;
            
            /* Select cells for each processor */
            forAll(multipleProcsPoints, iPoint)
            {
                gPointId = multipleProcsPoints[iPoint];
                lPointId = globalPointProcAddr[gPointId];
                const List<label>& lPointCells = cMesh_.pointCells(lPointId);
                gPointCells.resize(lPointCells.size());
                
                forAll(lPointCells, iCell)
                {
                    gPointCells[iCell] = localCellProcAddr[lPointCells[iCell]];
                    if (!cellProc_.found(gPointCells[iCell]))
                    {
                        cellProc_.insert(gPointCells[iCell], Pstream::myProcNo());
                    }
                }
                
                if (pointCells_.found(gPointId))
                {
                    pointCells_[gPointId].append(gPointCells);
                }
                else
                {
                    pointCells_.insert(gPointId,gPointCells);
                }
            }
            
            if (Pstream::master())
            {
                
                for (label proci = Pstream::firstSlave(); proci <= Pstream::lastSlave(); proci++)
                {
                    IPstream fromSlave(Pstream::commsTypes::scheduled, proci);
                    labelHashTable<List<label> > slaveCells;
                    labelHashTable<label> slaveCellProc;
                    
                    //gather pointCells_ table
                    fromSlave >> slaveCells;
                    
                    //copy data from slave to master array
                    forAllConstIter(labelHashTable<List<label> >, slaveCells, iter)
                    {
                        if (pointCells_.found(iter.key()))
                        {
                            pointCells_[iter.key()].append(iter());
                        }
                        else
                        {
                            pointCells_.insert(iter.key(), iter());
                        }
                    }
                    
                    //gather cellProc_ table
                    fromSlave >> slaveCellProc;
                    forAllConstIter(labelHashTable<label>, slaveCellProc, iter)
                    {
                        if (!cellProc_.found(iter.key()))
                        {
                            cellProc_.insert(iter.key(), iter());
                        }
                        else
                        {
                            if (iter() != proci)
                            {
                                Info << "Same cell located at two or more processors!" << endl;
                            }
                        }
                    }
                }
            }
            else
            {
                OPstream toMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
                toMaster << pointCells_;
                toMaster << cellProc_;
            }
        }
        
        
        //send accumulated data to slave processes
        if (Pstream::master())
        {
            for(label proci = Pstream::firstSlave(); proci <= Pstream::lastSlave(); proci++)
            {
                OPstream toSlave(Pstream::commsTypes::scheduled, proci);
                //toSlave << pointProcs_;
                toSlave << pointCells_;
                toSlave << cellProc_;
            }
        }
        else
        {
            //pointProcs_.clearStorage();
            pointCells_.clearStorage();
            cellProc_.clearStorage();
            
            IPstream fromMaster(Pstream::commsTypes::scheduled, Pstream::masterNo());
            //fromMaster >> pointProcs_;
            fromMaster >> pointCells_;
            fromMaster >> cellProc_;
        }
        
        Info << "Data redistributed" << endl;
        
        //loop over multiple processor vertices and find corresponding faces and patch id's
        
        labelHashTable<List<label> > faceIds;
        labelHashTable<label> faceProcPatch;
        labelHashTable<List<label> > faceNeigCells;
        labelHashTable<labelHashSet> neigProcFaces;
        List<DynamicList<label> > globalCorCellIds;
        labelHashTable<label> globalCorProcIds;
        //labelHashTable<List<label> > globalCellIds;
        
        {
            label gPointId = -1;
            label pointId  = -1;
            label faceId   = -1;
            label fPatchId = -1;
            //Pout << "multipleProcsPoints = " << multipleProcsPoints << endl;
            forAll(multipleProcsPoints, iPoint)
            {
                gPointId = multipleProcsPoints[iPoint];
                pointId  = globalPointProcAddr[gPointId];
                labelList locFaceIds = cMesh_.pointFaces()[pointId];
                forAll(locFaceIds, iFace)
                {
                    faceId = locFaceIds[iFace];
                    if (cMesh_.isInternalFace(faceId)) //do nothing within this implementation
                    {
                    }
                    else
                    {

                        bool isProcPatch = false;
                        fPatchId = cMesh_.boundaryMesh().whichPatch(faceId);
                        forAll(procPairs_, iPatch)
                        {
                            if (fPatchId == procPairs_[iPatch])
                            {
                                isProcPatch = true;
                                break;
                            }
                        }
                        if (isProcPatch)
                        {
                            if (faceIds.found(faceId))
                            {
                                faceIds[faceId].append(pointId);
                            }
                            else
                            {
                                faceIds.insert
                                (
                                    faceId,
                                    List<label>(1, pointId)
                                );
                                faceProcPatch.insert
                                (
                                    faceId,
                                    fPatchId
                                );
                            }
                        }
                    }
                }
            }
        }
        
        //store for each face - global cell id at neighbouring process
        {
            label lPointId = -1;
            label gPointId = -1;
            label gCellId  = -1;
            label procId   = -1;
            label patchId  = -1;
            forAllConstIter(labelHashTable<List<label> >, faceIds, iter)
            {
                List<label> pIds = iter();
                labelHashSet fCells;
                forAll (pIds, iPoint)
                {
                    lPointId = pIds[iPoint];
                    gPointId = localPointProcAddr[lPointId];
                    List<label> pCells = pointCells_[gPointId];
                    //store only uniqe cells and only from corner-neighbouring processes
                    forAll(pCells, iCell)
                    {
                        gCellId = pCells[iCell];
                        procId = cellProc_[gCellId];
                        //check that cell at "corner" process
                        //not self process or process connected through patch
                        bool cornerProcess = true;
                        if (Pstream::myProcNo() == procId)
                        {
                            cornerProcess = false;
                        }
                        patchId = cMesh_.boundaryMesh().whichPatch(iter.key());
                        if (neigProcs_[idProcPatchPairs_[patchId]] == procId)
                        {
                            cornerProcess = false;
                        }
                        
                        if (cornerProcess) // "corner" process
                        {
                            fCells.insert(gCellId); //make sure cell ids are unique
                            if (neigProcFaces.found(procId))
                            {
                                neigProcFaces[procId].insert(iter.key());
                            }
                            else
                            {
                                neigProcFaces.insert(procId, labelHashSet());
                                neigProcFaces[procId].insert(iter.key());
                            }
                        }
                    }
                }
                faceNeigCells.insert
                (
                    iter.key(),
                    fCells.toc()
                );
            }
        }
        
        //and create the list of cells ids, requested from neigb. process and addresing
        //table from these cells to current proc. patch faces
        {
            corAddr_.resize(neigProcFaces.toc().size());
            globalCorCellIds.resize(corAddr_.size());
            
            label iNeiProc = 0;
            label patchId  = -1;
            label patchNo  = -1;
            label faceId   = -1;
            label faceNo   = -1;
            label gCellId  = -1;
            labelHashTable<label> cellk;
            //Pout << "faceIds.toc() " << faceIds.toc() << endl;
            forAll(faceIds.toc(), iFaceId)
            {
                cellk.insert
                (
                    faceIds.toc()[iFaceId],
                    0
                );
            }
            
            forAllConstIter(labelHashTable<labelHashSet>, neigProcFaces, iter)
            {
                neigProcs_.append(iter.key());
                procPairs_.append(-1); //corner pair
                globalCorProcIds.insert(iter.key(),iNeiProc);
                
                List<label> procFaceIds = iter().toc();
                
                forAll(procFaceIds, iFace)
                {
                    faceId = procFaceIds[iFace];
                    patchId= cMesh_.boundaryMesh().whichPatch(faceId);
                    patchNo= idProcPatchPairs_[patchId];
                    faceNo = cMesh_.boundaryMesh()[patchId].whichFace(faceId);
                    
                    List<label> cellIds = faceNeigCells[faceId];
                    
                    forAll(cellIds, iCell)
                    {
                        gCellId = cellIds[iCell];
                        
                        if (cellProc_[gCellId] == iter.key())
                        {
                            globalCorCellIds[iNeiProc].append(gCellId);
                            corAddr_[iNeiProc].append
                            (
                                Triple<label>
                                ()
                            );
                            
                            corAddr_[iNeiProc].last()[0] = 
                                    patchNo;                      // # patch
                            corAddr_[iNeiProc].last()[1] = 
                                    faceNo;                       // # face at patch
                            corAddr_[iNeiProc].last()[2] = 
                                    cellk[faceId];                // # cell for face
                            cellk[faceId]++;
                        }
                    }
                    corEnd_[patchNo][faceNo] = corStart_[patchNo][faceNo] + cellk[faceId] - 1;
                }
                iNeiProc++; //next proc
            }
        }
        
        //redistribute global cell ids amongth other processors
        //and convert global addressing to local addressing
        forAll(neigProcs_, iProc)
        {
            if (procPairs_[iProc] < 0)
            {
                label procId = neigProcs_[iProc];
                label id = globalCorProcIds[procId];
                
                OPstream oStream (Pstream::commsTypes::scheduled, procId);
                oStream << globalCorCellIds[id];
            }
        }
        
        corCellIds_.resize(corAddr_.size());
        label iCorProc = 0;
        forAll(neigProcs_, iProc)
        {
            if (procPairs_[iProc] < 0)
            {
                label procId = neigProcs_[iProc];
                IPstream iStream (Pstream::commsTypes::scheduled, procId);
                List<label> neiGlobalCellIds (iStream);

                corProcIds_.insert
                (
                    procId,
                    iCorProc
                );

                corCellIds_[iCorProc].resize(neiGlobalCellIds.size());
                
                forAll(neiGlobalCellIds, iCell)
                {
                    corCellIds_[iCorProc][iCell] = 
                        globalCellProcAddr[neiGlobalCellIds[iCell]];
                }
                
                iCorProc++;
            }
        }
        
    }
    
};

//
//END-OF-FILE
//

