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
    Foam::fvsc::leastSquaresOpt::extendedFaceStencilFormVfValues
Description 
    Function for filling of array which contain which will contain
    the cells necessary for use by several processors.
Source file
    extendedFaceStencilScalarDer.C
    extendedFaceStencilScalarGradOpt.C
    leastSquaresStencilOpt.C
\*---------------------------------------------------------------------------*/

void formVfValues (const volScalarField& iF,List3<scalar>& procVfValues)
{
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
        //data size from processor, that is neighbouring through the patch
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
//            Pout << "corAddr for scalar: " << "corAddt["<< iCorProc << "] = " << corAddr_[iCorProc] << endl;

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

//    Info << "procVfValues = " << procVfValues << endl;
}


template<class FieldType>
void formVfValues (const GeometricField<FieldType, fvPatchField, volMesh>& iF,List<List3<scalar>>& procVfValues)
{
    //set values from this domain   
    label cellId = -1;
    label nComps = pTraits<FieldType>::nComponents;

    procVfValues.resize(nComps); 
    for(label compI = 0;compI <= nComps-1;compI++)
    {       
        forAll(procPairs_, patchI)
        {
            if (procPairs_[patchI] > -1)
            {
                procVfValues[compI].resize(procVfValues[compI].size()+1);
                procVfValues[compI][patchI].resize(procWf2_[patchI].size());

                forAll(procVfValues[compI][patchI], faceI)
                {
                    procVfValues[compI][patchI][faceI].resize(procWf2_[patchI][faceI].size());
//                    procVfValues[patchI][faceI] = 0.0; //make values zero
                    forAll(procVfValues[compI][patchI][faceI], cellI)
                    {
                        procVfValues[compI][patchI][faceI][cellI] =0.0; //make values zero              

//                        cellId = myProcPatchCells_[patchI][faceI][cellI];
                        procVfValues[compI][patchI][faceI][cellI] = iF.primitiveField()[cellId].component(compI);
                   
                    
                    }
                }
            }
        }
    }
 
    //Step 1. Send field data to neighbouring processors (non-blocking mode)

    PstreamBuffers pBuffers(Pstream::commsTypes::nonBlocking);

    forAll(procPairs_, procI)
    {
        label procId = neigProcs_[procI];
        List<DynamicList<scalar>> locVf;
        locVf.resize(nComps);

        for(label compI = 0;compI <= nComps-1;compI++)
        {   
            if (procPairs_[procI] > -1) //patch proc pair
            {        
                forAll(procVfValues[compI][procI], faceI)
                {
                    for(
                            label
                            cellI = 0;
                            cellI <= ownEnd_[procI][faceI];
                            cellI++
                       )
                    {
                    
                        locVf[compI].append(procVfValues[compI][procI][faceI][cellI]); 
                    }
                }
   
            }

            else //corner connected process
            {
                label cellId = -1;
                label addrId = corProcIds_[procId];
                forAll(corCellIds_[addrId], iCellId)
                {
                    cellId = corCellIds_[addrId][iCellId];
                    locVf[compI].append(iF.primitiveField()[cellId].component(compI));
                }
            }
        }
        UOPstream oProcStr(procId, pBuffers);
        oProcStr << locVf;
    }
    pBuffers.finishedSends();
    //Step 2. Recieve field data from neighbouring processors
    label iCorProc = 0;
    forAll(procPairs_, procI)
    {
//        data size from processor, that is neighbouring through the patch

        label procId = neigProcs_[procI];
        
        UIPstream iProcStr(procId, pBuffers);
         
        List<List<scalar>> locVf(iProcStr);  
        
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
                    for(label compI = 0;compI <= nComps-1;compI++)
                    procVfValues[compI][procI][iFace][iCell] =
                    locVf[compI][iVf];
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

//        Pout << "corAddr for template: " << "corAddt["<< iCorProc << "] = " << corAddr_[iCorProc] << endl; 
            forAll(addr, iVal)
            {
                patchNo = addr[iVal][0];
                faceNo  = addr[iVal][1];
                cellNo  = addr[iVal][2];
                offset = corStart_[patchNo][faceNo];          
         
                for(label compI = 0;compI <= nComps-1;compI++)
                procVfValues[compI][patchNo][faceNo][cellNo+offset] = locVf[compI][iVal];
            }
            iCorProc++;

        }
    }

};
