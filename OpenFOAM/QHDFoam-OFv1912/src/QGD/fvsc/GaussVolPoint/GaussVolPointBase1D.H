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
    Foam::fvsc::GaussVolPoint::GaussVolPointBase1D
SourceFile
    GaussVolPointBase1D.C
\*---------------------------------------------------------------------------*/

#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"
#include "tmp.H"

#ifndef GaussVolPointBase1D_H
#define GaussVolPointBase1D_H

namespace Foam
{

class fvMesh;

namespace fvsc
{

class GaussVolPointBase1D
{

    //-
    tmp<surfaceVectorField> nfRef_;

public:

    //-
    GaussVolPointBase1D(const fvMesh& mesh);
    
    //-
    ~GaussVolPointBase1D();
    
    //-
    void faceGrad(const volScalarField& f, surfaceVectorField& gradf);

    //-
    void faceGrad(const volVectorField& f, surfaceTensorField& gradf);

    //-
    void faceDiv(const volVectorField& f, surfaceScalarField& divf);
    
    //-
    void faceDiv(const volTensorField& f, surfaceVectorField& divf);

};

}

}

#endif

//
//END-OF-FILE
//

