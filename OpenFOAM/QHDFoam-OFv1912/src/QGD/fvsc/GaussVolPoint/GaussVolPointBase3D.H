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
    Foam::fvsc::GaussVolPoint::GaussVolPointBase3D
SourceFile
    GaussVolPointBase3D.C
\*---------------------------------------------------------------------------*/

#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "tmp.H"
#include "List.H"
#include "DynamicList.H"
#include "label.H"
#include "scalar.H"
#include "autoPtr.H"

#ifndef GaussVolPointBase3D_H
#define GaussVolPointBase3D_H

namespace Foam
{

class fvMesh;
class face;
class volPointInterpolation;

typedef List<face> faceList;

namespace fvsc
{

class GaussVolPointBase3D
{
    
    //-
    const volPointInterpolation& volPoint_;

    //+
    tmp<surfaceVectorField> nfRef_;

    //+
    List<List<label> > bgfid_;
    
    //+
    List<bool> processorPatch_;

    /* coefficients for quad faces */

    //+
    List<label> qf_;
    
    //+
    List<List<scalar> > aqx_;
    
    //+
    List<List<scalar> > aqy_;
    
    //+
    List<List<scalar> > aqz_;
    
    //+
    List<scalar> vq_;
    
    //+
    List<List<label> > bqf_;
    
    //+
    List<List<List<scalar> > > baqx_;

    //+
    List<List<List<scalar> > > baqy_;

    //+
    List<List<List<scalar> > > baqz_;
    
    //+
    List<List<scalar> > bvq_;
    
    /* coefficients for tri faces */
    //
    List<label> tf_;
    
    //
    List<List<scalar> > atx_;
    
    //
    List<List<scalar> > aty_;
    
    //
    List<List<scalar> > atz_;
    
    //
    List<scalar> vt_;
    
    //
    List<List<label> > btf_;
    
    //
    List<List<List<scalar> > > batx_;

    //
    List<List<List<scalar> > > baty_;

    //
    List<List<List<scalar> > > batz_;
    
    //-
    List<List<scalar> > bvt_;
    
    /* coefficients for other faces */
    
    //-
    List<List<scalar> > bmvON_;
    
    //+
    DynamicList<label> of_;
    
    //+
    List<List<label> > bof_;

protected:

    //-
    void calcGradfIF
    (
        const volScalarField& sf,
        const pointScalarField& pf,
        const faceList& faces,
        surfaceVectorField& gradf,
        const surfaceVectorField& dfdn
    );
    
    //-
    void calcGradfBF
    (
        const volScalarField& sf,
        const pointScalarField& pf,
        const faceList& faces,
        surfaceVectorField& gradf,
        const surfaceVectorField& dfdn
    );

    //-
    void calcGradfIF
    (
        const volVectorField& sf,
        const pointVectorField& pf,
        const faceList& faces,
        surfaceTensorField& gradf,
        const surfaceTensorField& dfdn
    );
    
    //-
    void calcGradfBF
    (
        const volVectorField& sf,
        const pointVectorField& pf,
        const faceList& faces,
        surfaceTensorField& gradf,
        const surfaceTensorField& dfdn
    );

    //-
    void calcDivfIF
    (
        const volVectorField& sf,
        const pointVectorField& pf,
        const faceList& faces,
        surfaceScalarField& divf,
        const surfaceScalarField& dfdn
    );

    //-
    void calcDivfBF
    (
        const volVectorField& sf,
        const pointVectorField& pf,
        const faceList& faces,
        surfaceScalarField& divf,
        const surfaceScalarField& dfdn
    );

    //-
    void calcDivfIF
    (
        const volTensorField& tf,
        const pointTensorField& pf,
        const faceList& faces,
        surfaceVectorField& divf,
        const surfaceVectorField& dfdn
    );

    //-
    void calcDivfBF
    (
        const volTensorField& sf,
        const pointTensorField& pf,
        const faceList& faces,
        surfaceVectorField& divf,
        const surfaceVectorField& dfdn
    );
    
    //- Calculate weights for triangles
    void triCalcWeights
    (
        const fvMesh& m
    );

    //- Calcualte weights for quads
    void quaCalcWeights
    (
        const fvMesh& m
    );

public:

    //-
    GaussVolPointBase3D(const fvMesh& mesh);
    
    //-
    ~GaussVolPointBase3D();
    
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

