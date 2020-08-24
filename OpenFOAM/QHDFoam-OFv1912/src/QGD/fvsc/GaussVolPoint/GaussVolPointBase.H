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
    Foam::fvsc::GaussVolPoint::GaussVolPointBase
SourceFile
    GaussVolPointBase.C
\*---------------------------------------------------------------------------*/

#ifndef GaussVolPointBase_H
#define GaussVolPointBase_H

#include "fvMesh.H"

#include "GaussVolPointBase1D.H"
#include "GaussVolPointBase2D.H"
#include "GaussVolPointBase3D.H"

namespace Foam
{

namespace fvsc
{

class GaussVolPointBase : 
    public GaussVolPointBase1D,
    public GaussVolPointBase2D,
    public GaussVolPointBase3D
{

public:

    //-
    GaussVolPointBase(const fvMesh& mesh);
    
    //-
    ~GaussVolPointBase();
    
    //-
    void faceGrad (const volScalarField& vF, surfaceVectorField& gradf);

    //-
    void faceGrad (const volVectorField& vF, surfaceTensorField& gradf);

    //-
    void faceDiv (const volVectorField& vVF, surfaceScalarField& divf);
    
    //-
    void faceDiv(const volTensorField& vTF, surfaceVectorField& divf);
};

}

}

#endif

//
//
//


