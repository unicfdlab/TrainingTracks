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
    Foam::fvsc::leastSquaresOpt::extendedFaceStencilScalarDer
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

void Foam::fvsc::leastSquaresOpt::faceScalarDer(const Field<scalar>& iF,const Field<scalar>& sF,int com, surfaceScalarField& rField)
{
    forAll(sF, facei)
    {
        rField[facei] = 0.0;
        forAll(GdfAll_[facei], i)
        {
            rField[facei] += wf2All_[facei][i]*GdfAll_[facei][i].component(com)*(iF[neighbourCells_[facei][i]] - sF[facei]);
        }
    }

};

void Foam::fvsc::leastSquaresOpt::faceScalarDer(const tmp<Field<scalar>>& tiF,const tmp<Field<scalar>>& tsF, int com, tmp<surfaceScalarField>& trField )
{
};


void Foam::fvsc::leastSquaresOpt::faceScalarDer(const List3<scalar>& procVfValues,const surfaceScalarField& sF,int derComp, surfaceScalarField& rField)
{   
    scalar gf = 0.0;
    forAll(procPairs_, patchI)
    {
        label procPatchId = procPairs_[patchI];
        if (procPatchId > -1)
        {
            fvsPatchScalarField& pgradf = rField.boundaryFieldRef()[procPatchId];
            const List2<scalar> & pvf = procVfValues[patchI];
            const List2<scalar> & pwf2= procWf2_[patchI];
            const List2<vector> & pgdf= procGdf_[patchI];

            forAll(pgradf, iFace)
            {
                gf = 0.0;
                forAll(procGdf_[patchI][iFace], i)
                {
                    gf += pwf2[iFace][i]*pgdf[iFace][i].component(derComp)*
                       (pvf[iFace][i] - sF.boundaryField()[procPatchId][iFace]);
                }
                rField[iFace] = gf;
            }
        }
    }
};

void Foam::fvsc::leastSquaresOpt::faceScalarDer(const tmp<List3<scalar>>& tprocVfValues,const tmp<surfaceScalarField>&, int derComp, tmp<surfaceScalarField>& trField )
{
};

