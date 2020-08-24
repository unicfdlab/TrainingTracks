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
    Foam::fvsc::GaussVolPoint::GaussVolPointBase2D
SourceFile
    GaussVolPointBase2D.C
\*---------------------------------------------------------------------------*/

#ifndef GaussVolPointBase2D_H
#define GaussVolPointBase2D_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"

namespace Foam
{

class fvMesh;

namespace fvsc
{

class GaussVolPointBase2D
{

protected:
    
    List<scalar> c1_;
    List<scalar> c2_;
    List<scalar> c3_;
    List<scalar> c4_;
//    List<scalar> cosa1_;
//    List<scalar> cosa2_;
//    List<scalar> sina1_;
//    List<scalar> sina2_;
//    List<scalar> den_;
//    List<vector> v42_;
//    List<vector> v13_;
    List<scalar> mv42_;
    List<scalar> mv13_;
    List<label>  ip3_;
    List<label>  ip1_;
    List<label>  ic4_;
    List<label>  ic2_;
    vector e1_;
    vector e2_;
    label ie1_, ie2_, ie3_;
    
    //- ordinary external boundaries
    //  or processor boundary
    List<label> ordinaryPatches_;
    List<bool>  processorPatch_;
    
//    List<List<scalar> > cosa1e_;
//    List<List<scalar> > cosa2e_;
//    List<List<scalar> > sina1e_;
//    List<List<scalar> > sina2e_;
    List<List<label> >  ip3e_;
    List<List<label> >  ip1e_;
    List<List<label> >  ic4e_;
    List<List<scalar> > c1e_;
    List<List<scalar> > c2e_;
    List<List<scalar> > c3e_;
    List<List<scalar> > c4e_;
    List<List<scalar> > mv42e_;
    List<List<scalar> > mv13e_;
//    List<List<scalar> > dene_;
    
public:
    
    //-
    GaussVolPointBase2D(const fvMesh& mesh);
    
    //-
    ~GaussVolPointBase2D();
    
    //-
    void faceGrad(const volScalarField& f, surfaceVectorField& gradf);
    
    //-
    void faceDiv (const volVectorField& f, surfaceScalarField& divf);
    
    //-
    void faceDiv (const volTensorField& f, surfaceVectorField& divf);
};

}

}

#endif

//
//
//


