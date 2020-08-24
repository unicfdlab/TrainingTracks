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
    Foam::fvsc::GaussVolPoint::GaussVolPoint3D
Description
    This is a method for approximating derivatives of tangents to a face (3D case). 
    They are further used in the calculation of the QGD terms.
\*---------------------------------------------------------------------------*/
#include "GaussVolPointBase3D.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"
#include "volPointInterpolation.H"
#include "processorFvPatch.H"
#include "processorFvPatchFields.H"

Foam::fvsc::GaussVolPointBase3D::GaussVolPointBase3D(const fvMesh& mesh)
:
    volPoint_
    (
        volPointInterpolation::New(mesh)
    ),
    nfRef_(mesh.thisDb().lookupObject<surfaceVectorField>("nf")),
    bgfid_(mesh.boundary().size()),
    processorPatch_(mesh.boundary().size(), false),
    qf_(0),
    aqx_(0),
    aqy_(0),
    aqz_(0),
    vq_(0),
    bqf_(mesh.boundary().size()),
    baqx_(mesh.boundary().size()),
    baqy_(mesh.boundary().size()),
    baqz_(mesh.boundary().size()),
    bvq_(mesh.boundary().size()),
    tf_(0),
    atx_(0),
    aty_(0),
    atz_(0),
    vt_(0),
    btf_(mesh.boundary().size()),
    batx_(mesh.boundary().size()),
    baty_(mesh.boundary().size()),
    batz_(mesh.boundary().size()),
    bvt_(mesh.boundary().size()),
    bmvON_(mesh.boundary().size()),
    of_(0),
    bof_(mesh.boundary().size())
{
    forAll(bgfid_, iPatch)
    {
        if (isA<processorFvPatch>(mesh.boundary()[iPatch]))
        {
            processorPatch_[iPatch] = true;
        }
        const fvPatch& fvp = mesh.boundary()[iPatch];
        bgfid_[iPatch].resize(fvp.size());
        forAll(fvp, i)
        {
            bgfid_[iPatch][i] = mesh.boundary()[iPatch].start() + i;
        }
    }
    
    //sort quad, tri faces and other
    const faceList& faces = mesh.faces();
    
    forAll(faces, i)
    {
        if (mesh.isInternalFace(i))
        {
            if (faces[i].size() == 3)
            {
                tf_.append(i);
            }
            else if (faces[i].size() == 4)
            {
                qf_.append(i);
            }
            else
            {
                of_.append(i);
            }
        }
    }
    forAll(bgfid_, iPatch)
    {
        const fvPatch& fvp = mesh.boundary()[iPatch];
        forAll(fvp, i)
        {
            if (faces[bgfid_[iPatch][i]].size() == 3)
            {
                btf_[iPatch].append(i);
            }
            else if (faces[bgfid_[iPatch][i]].size() == 4)
            {
                bqf_[iPatch].append(i);
            }
            else
            {
                bof_[iPatch].append(i);
            }
        }
    }
    
    forAll(bgfid_, iPatch)
    {
        vectorField vO(mesh.boundary()[iPatch].size(), vector::zero);
        vectorField vN(mesh.boundary()[iPatch].size(), vector::zero);
        
        vO = mesh.boundary()[iPatch].Cn();
        if (processorPatch_[iPatch])
        {
            vN = refCast<const processorFvPatch>(mesh.boundary()[iPatch]).
                    procPolyPatch().neighbFaceCellCentres();
        }
        else
        {
            vN = vO + 2.0*
            (
                mesh.boundary()[iPatch].Cf()
                -
                vO
            );
        }
        
        bmvON_[iPatch] = mag
        (
            vO - vN
        );
    }

    triCalcWeights(mesh);
    
    quaCalcWeights(mesh);
};

void Foam::fvsc::GaussVolPointBase3D::triCalcWeights(const fvMesh& mesh)
{
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    label facei = -1;
    atx_.resize(tf_.size());
    aty_.resize(tf_.size());
    atz_.resize(tf_.size());
    vt_.resize(tf_.size());
    label own = -1;
    label nei = -1;
    label p1 = -1, p2 = -1, p3 = -1;
    const scalar OneBySix = (1.0 / 6.0);
    
    forAll(tf_, i)
    {
        facei = tf_[i];
        own   = mesh.owner()[facei];
        nei   = mesh.neighbour()[facei];
        p1    = faces[facei][0];
        p2    = faces[facei][1];
        p3    = faces[facei][2];
        //p4 - is nei
        //p5 - is own

        vt_[i] =
        (
            (points[p2] - points[p1]) ^ (points[p3] - points[p1])
        ) & (mesh.C()[own] - mesh.C()[nei]);
        vt_[i] *= OneBySix;
        
        /* Coefficients for X */
        atx_[i].resize(5);
        atx_[i][0] = OneBySix*((mesh.C()[own].z() - mesh.C()[nei].z())*(points[p2].y() - points[p3].y()) + 
            (mesh.C()[nei].y() - mesh.C()[own].y())*(points[p2].z() - points[p3].z()));
        atx_[i][1] = OneBySix*((mesh.C()[nei].y() - mesh.C()[own].y())*(points[p3].z() - points[p1].z()) + 
            (mesh.C()[own].z() - mesh.C()[nei].z())*(points[p3].y() - points[p1].y()));
        atx_[i][2] = OneBySix*((mesh.C()[nei].y() - mesh.C()[own].y())*(points[p1].z() - points[p2].z()) + 
            (mesh.C()[own].z() - mesh.C()[nei].z())*(points[p1].y() - points[p2].y()));
        atx_[i][3] = OneBySix*(points[p1].z()*(points[p2].y() - points[p3].y()) + 
                               points[p2].z()*(points[p3].y() - points[p1].y()) + 
                               points[p3].z()*(points[p1].y() - points[p2].y()));
        atx_[i][4] = -atx_[i][3];
        
        /* Coefficients for Y */
        aty_[i].resize(5);
        aty_[i][0] = OneBySix*((mesh.C()[own].x() - mesh.C()[nei].x())*(points[p2].z() - points[p3].z()) + 
            (mesh.C()[nei].z() - mesh.C()[own].z())*(points[p2].x() - points[p3].x()));
        aty_[i][1] = OneBySix*((mesh.C()[nei].z() - mesh.C()[own].z())*(points[p3].x() - points[p1].x()) + 
            (mesh.C()[own].x() - mesh.C()[nei].x())*(points[p3].z() - points[p1].z()));
        aty_[i][2] = OneBySix*((mesh.C()[nei].z() - mesh.C()[own].z())*(points[p1].x() - points[p2].x()) + 
            (mesh.C()[own].x() - mesh.C()[nei].x())*(points[p1].z() - points[p2].z()));
        aty_[i][3] = OneBySix*(points[p1].x()*(points[p2].z() - points[p3].z()) + 
                               points[p2].x()*(points[p3].z() - points[p1].z()) + 
                               points[p3].x()*(points[p1].z() - points[p2].z()));
        aty_[i][4] = -aty_[i][3];
        
        /* Coefficients for Z */
        atz_[i].resize(5);
        atz_[i][0] = OneBySix*((mesh.C()[own].y() - mesh.C()[nei].y())*(points[p2].x() - points[p3].x()) + 
            (mesh.C()[nei].x() - mesh.C()[own].x())*(points[p2].y() - points[p3].y()));
        atz_[i][1] = OneBySix*((mesh.C()[nei].x() - mesh.C()[own].x())*(points[p3].y() - points[p1].y()) + 
            (mesh.C()[own].y() - mesh.C()[nei].y())*(points[p3].x() - points[p1].x()));
        atz_[i][2] = OneBySix*((mesh.C()[nei].x() - mesh.C()[own].x())*(points[p1].y() - points[p2].y()) + 
            (mesh.C()[own].y() - mesh.C()[nei].y())*(points[p1].x() - points[p2].x()));
        atz_[i][3] = OneBySix*(points[p1].y()*(points[p2].x() - points[p3].x()) + 
                               points[p2].y()*(points[p3].x() - points[p1].x()) + 
                               points[p3].y()*(points[p1].x() - points[p2].x()));
        atz_[i][4] = -atz_[i][3];
    }
    
    forAll(btf_, iPatch)
    {
        batx_[iPatch].resize(btf_[iPatch].size());
        baty_[iPatch].resize(btf_[iPatch].size());
        batz_[iPatch].resize(btf_[iPatch].size());
        bvt_[iPatch].resize(btf_[iPatch].size());
        
        vectorField v5(mesh.boundary()[iPatch].size(), vector::zero);
        vectorField v4(mesh.boundary()[iPatch].size(), vector::zero);
        
        v5 = mesh.boundary()[iPatch].Cn();
        if (processorPatch_[iPatch])
        {
            v4 = refCast<const processorFvPatch>(mesh.boundary()[iPatch]).
                    procPolyPatch().neighbFaceCellCentres();
        }
        else
        {
            v4 = v5 + 2.0*
            (
                mesh.boundary()[iPatch].Cf()
                -
                v5
            );
        }
        
        label gFaceId = -1;
        
        forAll(btf_[iPatch], k)
        {
            facei = btf_[iPatch][k];
            gFaceId = bgfid_[iPatch][facei];
            
            p1  = faces[gFaceId][0];
            p2  = faces[gFaceId][1];
            p3  = faces[gFaceId][2];
            
            //p4 - is nei and stored in v4
            //p5 - is own and stored in v5
            
            bvt_[iPatch][k] =
                ((points[p2] - points[p1]) ^ (points[p3] - points[p1]))
                & (v5[facei] - v4[facei]);
            //bvt_[iPatch][k] *= TwoBySix;
            bvt_[iPatch][k] *= OneBySix;
            
            /* Coefficients for X */
            batx_[iPatch][k].resize(5);
            batx_[iPatch][k][0] = OneBySix*((v5[facei].z() - v4[facei].z())*(points[p2].y() - points[p3].y()) + 
                (v4[facei].y() - v5[facei].y())*(points[p2].z() - points[p3].z()));
            batx_[iPatch][k][1] = OneBySix*((v4[facei].y() - v5[facei].y())*(points[p3].z() - points[p1].z()) + 
                (v5[facei].z() - v4[facei].z())*(points[p3].y() - points[p1].y()));
            batx_[iPatch][k][2] = OneBySix*((v4[facei].y() - v5[facei].y())*(points[p1].z() - points[p2].z()) + 
                (v5[facei].z() - v4[facei].z())*(points[p1].y() - points[p2].y()));
            batx_[iPatch][k][3] = OneBySix*(points[p1].z()*(points[p2].y() - points[p3].y()) + 
                                            points[p2].z()*(points[p3].y() - points[p1].y()) + 
                                            points[p3].z()*(points[p1].y() - points[p2].y()));
            batx_[iPatch][k][4] = -batx_[iPatch][k][3];
            
            /* Coefficients for Y */
            baty_[iPatch][k].resize(5);
            baty_[iPatch][k][0] = OneBySix*((v5[facei].x() - v4[facei].x())*(points[p2].z() - points[p3].z()) + 
                (v4[facei].z() - v5[facei].z())*(points[p2].x() - points[p3].x()));
            baty_[iPatch][k][1] = OneBySix*((v4[facei].z() - v5[facei].z())*(points[p3].x() - points[p1].x()) + 
                (v5[facei].x() - v4[facei].x())*(points[p3].z() - points[p1].z()));
            baty_[iPatch][k][2] = OneBySix*((v4[facei].z() - v5[facei].z())*(points[p1].x() - points[p2].x()) + 
                (v5[facei].x() - v4[facei].x())*(points[p1].z() - points[p2].z()));
            baty_[iPatch][k][3] = OneBySix*(points[p1].x()*(points[p2].z() - points[p3].z()) + 
                                            points[p2].x()*(points[p3].z() - points[p1].z()) + 
                                            points[p3].x()*(points[p1].z() - points[p2].z()));
            baty_[iPatch][k][4] = -baty_[iPatch][k][3];
            
            /* Coefficients for Z */
            batz_[iPatch][k].resize(5);
            batz_[iPatch][k][0] = OneBySix*((v5[facei].y() - v4[facei].y())*(points[p2].x() - points[p3].x()) + 
                (v4[facei].x() - v5[facei].x())*(points[p2].y() - points[p3].y()));
            batz_[iPatch][k][1] = OneBySix*((v4[facei].x() - v5[facei].x())*(points[p3].y() - points[p1].y()) + 
                (v5[facei].y() - v4[facei].y())*(points[p3].x() - points[p1].x()));
            batz_[iPatch][k][2] = OneBySix*((v4[facei].x() - v5[facei].x())*(points[p1].y() - points[p2].y()) + 
                (v5[facei].y() - v4[facei].y())*(points[p1].x() - points[p2].x()));
            batz_[iPatch][k][3] = OneBySix*(points[p1].y()*(points[p2].x() - points[p3].x()) + 
                                            points[p2].y()*(points[p3].x() - points[p1].x()) + 
                                            points[p3].y()*(points[p1].x() - points[p2].x()));
            batz_[iPatch][k][4] = -batz_[iPatch][k][3];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::quaCalcWeights(const fvMesh& mesh)
{
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    const scalar OneBySix = (1.0 / 6.0);
    label facei = -1;
    aqx_.resize(qf_.size());
    aqy_.resize(qf_.size());
    aqz_.resize(qf_.size());
    vq_.resize(qf_.size());
    label own = -1;
    label nei = -1;
    label p4 = -1, p1 = -1, p2 = -1, p3 = -1;

    forAll(qf_, i)
    {
        facei = qf_[i];
        own = mesh.owner()[facei];
        nei = mesh.neighbour()[facei];
        p1  = faces[facei][0];
        p2  = faces[facei][1];
        p3  = faces[facei][2];
        p4  = faces[facei][3];
        //p5 - is nei
        //p6 - is own
        
        vq_[i] = (points[p3] - points[p1]) & 
        (
            (points[p4] - points[p2]) ^ (mesh.C()[own] - mesh.C()[nei])
        );
        vq_[i] *= OneBySix;
        
        /* coefficient for X */
        aqx_[i].resize(6);
        aqx_[i][0] = OneBySix*((mesh.C()[nei].y() - mesh.C()[own].y())*(points[p2].z() - points[p4].z())
            - (mesh.C()[nei].z() - mesh.C()[own].z())*(points[p2].y() - points[p4].y()));
        aqx_[i][1] = OneBySix*((mesh.C()[nei].y() - mesh.C()[own].y())*(points[p3].z() - points[p1].z())
            - (mesh.C()[nei].z() - mesh.C()[own].z())*(points[p3].y() - points[p1].y()));
        aqx_[i][5] = OneBySix*((points[p1].y() - points[p3].y())*(points[p2].z() - points[p4].z())
            - (points[p1].z() - points[p3].z())*(points[p2].y() - points[p4].y()));
        
        aqx_[i][2] = -aqx_[i][0];
        aqx_[i][3] = -aqx_[i][1];
        aqx_[i][4] = -aqx_[i][5];
        
        /* coefficient for Y */
        aqy_[i].resize(6);
        aqy_[i][0] = OneBySix*((mesh.C()[nei].z() - mesh.C()[own].z())*(points[p2].x() - points[p4].x())
            - (mesh.C()[nei].x() - mesh.C()[own].x())*(points[p2].z() - points[p4].z()));
        aqy_[i][1] = OneBySix*((mesh.C()[nei].z() - mesh.C()[own].z())*(points[p3].x() - points[p1].x())
            - (mesh.C()[nei].x() - mesh.C()[own].x())*(points[p3].z() - points[p1].z()));
        aqy_[i][5] = OneBySix*((points[p1].z() - points[p3].z())*(points[p2].x() - points[p4].x())
            - (points[p1].x() - points[p3].x())*(points[p2].z() - points[p4].z()));
        
        aqy_[i][2] = -aqy_[i][0];
        aqy_[i][3] = -aqy_[i][1];
        aqy_[i][4] = -aqy_[i][5];
        
        /* coefficient for Z */
        aqz_[i].resize(6);
        aqz_[i][0] = OneBySix*((mesh.C()[nei].x() - mesh.C()[own].x())*(points[p2].y() - points[p4].y())
            - (mesh.C()[nei].y() - mesh.C()[own].y())*(points[p2].x() - points[p4].x()));
        aqz_[i][1] = OneBySix*((mesh.C()[nei].x() - mesh.C()[own].x())*(points[p3].y() - points[p1].y())
            - (mesh.C()[nei].y() - mesh.C()[own].y())*(points[p3].x() - points[p1].x()));
        aqz_[i][5] = OneBySix*((points[p1].x() - points[p3].x())*(points[p2].y() - points[p4].y())
            - (points[p1].y() - points[p3].y())*(points[p2].x() - points[p4].x()));
        
        aqz_[i][2] = -aqz_[i][0];
        aqz_[i][3] = -aqz_[i][1];
        aqz_[i][4] = -aqz_[i][5];
    }
    forAll(bqf_, iPatch)
    {
        baqx_[iPatch].resize(bqf_[iPatch].size());
        baqy_[iPatch].resize(bqf_[iPatch].size());
        baqz_[iPatch].resize(bqf_[iPatch].size());
        bvq_[iPatch].resize(bqf_[iPatch].size());
        
        vectorField v6(mesh.boundary()[iPatch].size(), vector::zero);
        vectorField v5(mesh.boundary()[iPatch].size(), vector::zero);
        
        v6 = mesh.boundary()[iPatch].Cn();
        if (processorPatch_[iPatch])
        {
            v5 = refCast<const processorFvPatch>(mesh.boundary()[iPatch]).
                    procPolyPatch().neighbFaceCellCentres();
        }
        else
        {
            v5 = v6 + 2.0*
            (
                mesh.boundary()[iPatch].Cf()
                -
                v6
            );
        }
        
        label gFaceId = -1;
        forAll(bqf_[iPatch], k)
        {
            facei = bqf_[iPatch][k];
            gFaceId = bgfid_[iPatch][facei];
            
            p1  = faces[gFaceId][0];
            p2  = faces[gFaceId][1];
            p3  = faces[gFaceId][2];
            p4  = faces[gFaceId][3];
            //p5 - is nei and stored in v5
            //p6 - is own and stored in v6
            
            bvq_[iPatch][k] = (points[p3] - points[p1]) & 
            (
                (points[p4] - points[p2]) ^ (v6[facei] - v5[facei])
            );
            bvq_[iPatch][k] *= OneBySix;
            
            /* coefficient for X */
            baqx_[iPatch][k].resize(6);
            baqx_[iPatch][k][0] = OneBySix*((v5[facei].y() - v6[facei].y())*(points[p2].z() - points[p4].z())
                - (v5[facei].z() - v6[facei].z())*(points[p2].y() - points[p4].y()));
            baqx_[iPatch][k][1] = OneBySix*((v5[facei].y() - v6[facei].y())*(points[p3].z() - points[p1].z())
                - (v5[facei].z() - v6[facei].z())*(points[p3].y() - points[p1].y()));
            baqx_[iPatch][k][5] = OneBySix*((points[p1].y() - points[p3].y())*(points[p2].z() - points[p4].z())
                - (points[p1].z() - points[p3].z())*(points[p2].y() - points[p4].y()));
            
            baqx_[iPatch][k][2] = -baqx_[iPatch][k][0];
            baqx_[iPatch][k][3] = -baqx_[iPatch][k][1];
            baqx_[iPatch][k][4] = -baqx_[iPatch][k][5];

            /* coefficient for Y */
            baqy_[iPatch][k].resize(6);
            baqy_[iPatch][k][0] = OneBySix*((v5[facei].z() - v6[facei].z())*(points[p2].x() - points[p4].x())
                - (v5[facei].x() - v6[facei].x())*(points[p2].z() - points[p4].z()));
            baqy_[iPatch][k][1] = OneBySix*((v5[facei].z() - v6[facei].z())*(points[p3].x() - points[p1].x())
                - (v5[facei].x() - v6[facei].x())*(points[p3].z() - points[p1].z()));
            baqy_[iPatch][k][5] = OneBySix*((points[p1].z() - points[p3].z())*(points[p2].x() - points[p4].x())
                - (points[p1].x() - points[p3].x())*(points[p2].z() - points[p4].z()));
            
            baqy_[iPatch][k][2] = -baqy_[iPatch][k][0];
            baqy_[iPatch][k][3] = -baqy_[iPatch][k][1];
            baqy_[iPatch][k][4] = -baqy_[iPatch][k][5];
            
            /* coefficient for Z */
            baqz_[iPatch][k].resize(6);
            baqz_[iPatch][k][0] = OneBySix*((v5[facei].x() - v6[facei].x())*(points[p2].y() - points[p4].y())
                - (v5[facei].y() - v6[facei].y())*(points[p2].x() - points[p4].x()));
            baqz_[iPatch][k][1] = OneBySix*((v5[facei].x() - v6[facei].x())*(points[p3].y() - points[p1].y())
                - (v5[facei].y() - v6[facei].y())*(points[p3].x() - points[p1].x()));
            baqz_[iPatch][k][5] = OneBySix*((points[p1].x() - points[p3].x())*(points[p2].y() - points[p4].y())
                - (points[p1].y() - points[p3].y())*(points[p2].x() - points[p4].x()));
            
            baqz_[iPatch][k][2] = -baqz_[iPatch][k][0];
            baqz_[iPatch][k][3] = -baqz_[iPatch][k][1];
            baqz_[iPatch][k][4] = -baqz_[iPatch][k][5];
        }
    }
}

Foam::fvsc::GaussVolPointBase3D::~GaussVolPointBase3D()
{
}

#define VEC_CMPT(V,CMPT)\
    V[CMPT]

#define SCA_CMPT(V,CMPT)\
    V

#define dfdxif(vf,pf,dfdxfield,fi,vi,ai,icmpt,ocmpt,iop,oop)            \
{                                                                       \
    label celll = -1;                                                   \
    label facei = -1;                                                   \
    const label iown  = (ai.size() > 0) ? ai[0].size() - 1 : 0;         \
    const label inei  = (ai.size() > 0) ? iown - 1 : 0;                 \
    scalar dfdxface = 0.0;                                              \
    forAll(fi, i)                                                       \
    {                                                                   \
        facei = fi[i];                                                  \
        celll = vf.mesh().neighbour()[facei];                           \
        dfdxface =                                                      \
            iop(vf.primitiveField()[celll],icmpt) * ai[i][inei];        \
        celll = vf.mesh().owner()[facei];                               \
        dfdxface +=                                                     \
            iop(vf.primitiveField()[celll],icmpt) * ai[i][iown];        \
                                                                        \
        forAll(faces[facei], k)                                         \
        {                                                               \
            dfdxface +=                                                 \
                iop(pf[faces[facei][k]],icmpt) * ai[i][k];              \
        }                                                               \
        oop(dfdxfield.primitiveFieldRef()[facei],ocmpt)                 \
            += (dfdxface / vi[i]);                                      \
    }                                                                   \
}

#define dfdxbf(vf,pf,patchi,dfdxfield,bfi,bvfi,bai,icmpt,ocmpt,iop,oop) \
{                                                                       \
    label ifacei = -1;                                                  \
    const label iown  = (bai[patchi].size() > 0) ? bai[patchi][0].size() - 1 : 0; \
    const label inei  = (bai[patchi].size() > 0) ? iown - 1 : 0;                  \
    label gFaceId = -1;                                                 \
    scalar dfdxface = 0.0;                                              \
    forAll(bfi[patchi], i)                                              \
    {                                                                   \
        ifacei = bfi[patchi][i];                                        \
        gFaceId = bgfid_[patchi][ifacei];                               \
        dfdxface =                                                      \
            iop(psin[patchi][ifacei],icmpt) * bai[patchi][i][inei];     \
        dfdxface +=                                                     \
            iop(psio[ifacei],icmpt) * bai[patchi][i][iown];             \
        forAll(faces[gFaceId], k)                                       \
        {                                                               \
            dfdxface+=                                                  \
                    iop(pf[faces[gFaceId][k]],icmpt) *                  \
                    bai[patchi][i][k];                                  \
        }                                                               \
        oop(dfdxfield.boundaryFieldRef()[patchi][ifacei],ocmpt) +=      \
        dfdxface / bvfi[patchi][i];                                     \
    }                                                                   \
}



void Foam::fvsc::GaussVolPointBase3D::calcDivfIF
(
    const volVectorField& vf,
    const pointVectorField& pf,
    const faceList& faces,
    surfaceScalarField& divf,
    const surfaceScalarField& dfdn
)
{
    //calculate at quad faces
    dfdxif(vf,pf,divf,qf_,vq_,aqx_,0,0,VEC_CMPT,SCA_CMPT) //X
    dfdxif(vf,pf,divf,qf_,vq_,aqy_,1,0,VEC_CMPT,SCA_CMPT) //Y
    dfdxif(vf,pf,divf,qf_,vq_,aqz_,2,0,VEC_CMPT,SCA_CMPT) //Z
    
    //calculate at triangular faces
    dfdxif(vf,pf,divf,tf_,vt_,atx_,0,0,VEC_CMPT,SCA_CMPT) //X
    dfdxif(vf,pf,divf,tf_,vt_,aty_,1,0,VEC_CMPT,SCA_CMPT) //Y
    dfdxif(vf,pf,divf,tf_,vt_,atz_,2,0,VEC_CMPT,SCA_CMPT) //Z

    //other faces
    {
        label facei = -1;
        forAll(of_, i)
        {
            facei = of_[i];
            divf.primitiveFieldRef()[facei] = 
                dfdn.primitiveField()[facei];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcDivfBF
(
    const volVectorField& vf,
    const pointVectorField& pf,
    const faceList& faces,
    surfaceScalarField& divf,
    const surfaceScalarField& dfdn
)
{
    List<List<vector> > psin (vf.boundaryField().size());
    forAll(psin, iPatch)
    {
        if (processorPatch_[iPatch])
        {
            psin[iPatch] = refCast<const processorFvPatchField<vector> >
                (vf.boundaryField()[iPatch]).patchNeighbourField();
        }
        else
        {
            psin[iPatch] = vf.boundaryField()[iPatch] + 
                vf.boundaryField()[iPatch].snGrad()
                *
                bmvON_[iPatch]*0.5;
        }
        
        vectorField psio (vf.boundaryField()[iPatch].patchInternalField());
        
        // quad faces
        dfdxbf(vf,pf,iPatch,divf,bqf_,bvq_,baqx_,0,0,VEC_CMPT,SCA_CMPT) //X
        dfdxbf(vf,pf,iPatch,divf,bqf_,bvq_,baqy_,1,0,VEC_CMPT,SCA_CMPT) //Y
        dfdxbf(vf,pf,iPatch,divf,bqf_,bvq_,baqz_,2,0,VEC_CMPT,SCA_CMPT) //Z

        // tri faces
        dfdxbf(vf,pf,iPatch,divf,btf_,bvt_,batx_,0,0,VEC_CMPT,SCA_CMPT) //X
        dfdxbf(vf,pf,iPatch,divf,btf_,bvt_,baty_,1,0,VEC_CMPT,SCA_CMPT) //Y
        dfdxbf(vf,pf,iPatch,divf,btf_,bvt_,batz_,2,0,VEC_CMPT,SCA_CMPT) //Z

        //for other faces - apply surface normal derivative
        {
            label facei = -1;
            forAll(bof_[iPatch], i)
            {
                facei = bof_[iPatch][i];
                divf.boundaryFieldRef()[iPatch][facei] = 
                    dfdn.boundaryField()[iPatch][facei];
            }
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcDivfIF
(
    const volTensorField& tf,
    const pointTensorField& pf,
    const faceList& faces,
    surfaceVectorField& divf,
    const surfaceVectorField& dfdn
)
{
    //calculate at quandrangle faces
    //X
    dfdxif(tf,pf,divf,qf_,vq_,aqx_,0,0,VEC_CMPT,VEC_CMPT) // dT_xx / dx
    dfdxif(tf,pf,divf,qf_,vq_,aqy_,3,0,VEC_CMPT,VEC_CMPT) // dT_yx / dy
    dfdxif(tf,pf,divf,qf_,vq_,aqz_,6,0,VEC_CMPT,VEC_CMPT) // dT_zx / dz
    //Y
    dfdxif(tf,pf,divf,qf_,vq_,aqx_,1,1,VEC_CMPT,VEC_CMPT) // dT_xy / dx
    dfdxif(tf,pf,divf,qf_,vq_,aqy_,4,1,VEC_CMPT,VEC_CMPT) // dT_yy / dy
    dfdxif(tf,pf,divf,qf_,vq_,aqz_,7,1,VEC_CMPT,VEC_CMPT) // dT_zy / dz
    //Z
    dfdxif(tf,pf,divf,qf_,vq_,aqx_,2,2,VEC_CMPT,VEC_CMPT) // dT_xz / dx
    dfdxif(tf,pf,divf,qf_,vq_,aqy_,5,2,VEC_CMPT,VEC_CMPT) // dT_yz / dy
    dfdxif(tf,pf,divf,qf_,vq_,aqz_,8,2,VEC_CMPT,VEC_CMPT) // dT_zz / dz

    //calculate at triangular faces
    //X
    dfdxif(tf,pf,divf,tf_,vt_,atx_,0,0,VEC_CMPT,VEC_CMPT) // dT_xx / dx
    dfdxif(tf,pf,divf,tf_,vt_,aty_,3,0,VEC_CMPT,VEC_CMPT) // dT_yx / dy
    dfdxif(tf,pf,divf,tf_,vt_,atz_,6,0,VEC_CMPT,VEC_CMPT) // dT_zx / dz
    //Y
    dfdxif(tf,pf,divf,tf_,vt_,atx_,1,1,VEC_CMPT,VEC_CMPT) // dT_xy / dx
    dfdxif(tf,pf,divf,tf_,vt_,aty_,4,1,VEC_CMPT,VEC_CMPT) // dT_yy / dy
    dfdxif(tf,pf,divf,tf_,vt_,atz_,7,1,VEC_CMPT,VEC_CMPT) // dT_zy / dz
    //Z
    dfdxif(tf,pf,divf,tf_,vt_,atx_,2,2,VEC_CMPT,VEC_CMPT) // dT_xz / dx
    dfdxif(tf,pf,divf,tf_,vt_,aty_,5,2,VEC_CMPT,VEC_CMPT) // dT_yz / dy
    dfdxif(tf,pf,divf,tf_,vt_,atz_,8,2,VEC_CMPT,VEC_CMPT) // dT_zz / dz

    //
    //other faces
    {
        label facei = -1;
        forAll(of_, i)
        {
            facei = of_[i];
            divf.primitiveFieldRef()[facei] = 
                dfdn.primitiveField()[facei];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcDivfBF
(
    const volTensorField& tf,
    const pointTensorField& pf,
    const faceList& faces,
    surfaceVectorField& divf,
    const surfaceVectorField& dfdn
)
{
    List<List<tensor> > psin (tf.boundaryField().size());
    forAll(psin, iPatch)
    {
        if (processorPatch_[iPatch])
        {
            psin[iPatch] = refCast<const processorFvPatchField<tensor> >
                (tf.boundaryField()[iPatch]).patchNeighbourField();
        }
        else
        {
            psin[iPatch] = tf.boundaryField()[iPatch] + 
                tf.boundaryField()[iPatch].snGrad()
                *
                bmvON_[iPatch]*0.5;
        }
        
        tensorField psio (tf.boundaryField()[iPatch].patchInternalField());
        
        // quad faces
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqx_,0,0,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqy_,3,0,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqz_,6,0,VEC_CMPT,VEC_CMPT) //d/dz
        
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqx_,1,1,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqy_,4,1,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqz_,7,1,VEC_CMPT,VEC_CMPT) //d/dz
        
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqx_,2,2,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqy_,5,2,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,bqf_,bvq_,baqz_,8,2,VEC_CMPT,VEC_CMPT) //d/dz
        
        // tri faces
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,batx_,0,0,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,baty_,3,0,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,batz_,6,0,VEC_CMPT,VEC_CMPT) //d/dz
        
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,batx_,1,1,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,baty_,4,1,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,batz_,7,1,VEC_CMPT,VEC_CMPT) //d/dz
        
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,batx_,2,2,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,baty_,5,2,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,btf_,bvt_,batz_,8,2,VEC_CMPT,VEC_CMPT) //d/dz
        
        //for other faces - apply surface normal derivative
        {
            label facei = -1;
            forAll(bof_[iPatch], i)
            {
                facei = bof_[iPatch][i];
                divf.boundaryFieldRef()[iPatch][facei] = 
                    dfdn.boundaryField()[iPatch][facei];
            }
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcGradfIF
(
    const volScalarField& sf,
    const pointScalarField& pf,
    const faceList& faces,
    surfaceVectorField& gradf,
    const surfaceVectorField& dfdn
)
{
    //quad faces
    dfdxif(sf,pf,gradf,qf_,vq_,aqx_,0,0,SCA_CMPT,VEC_CMPT) //X
    dfdxif(sf,pf,gradf,qf_,vq_,aqy_,0,1,SCA_CMPT,VEC_CMPT) //Y
    dfdxif(sf,pf,gradf,qf_,vq_,aqz_,0,2,SCA_CMPT,VEC_CMPT) //Z
    
    //triangular faces
    dfdxif(sf,pf,gradf,tf_,vt_,atx_,0,0,SCA_CMPT,VEC_CMPT) //X
    dfdxif(sf,pf,gradf,tf_,vt_,aty_,0,1,SCA_CMPT,VEC_CMPT) //Y
    dfdxif(sf,pf,gradf,tf_,vt_,atz_,0,2,SCA_CMPT,VEC_CMPT) //Z

    //other faces
    {
        label facei = -1;
        forAll(of_, i)
        {
            facei = of_[i];
            gradf.primitiveFieldRef()[facei] = 
                dfdn.primitiveField()[facei];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcGradfBF
(
    const volScalarField& sf,
    const pointScalarField& pf,
    const faceList& faces,
    surfaceVectorField& gradf,
    const surfaceVectorField& dfdn
)
{
    List<List<scalar> > psin (sf.boundaryField().size());
    forAll(psin, iPatch)
    {
        if (processorPatch_[iPatch])
        {
            psin[iPatch] = refCast<const processorFvPatchField<scalar> >
                (sf.boundaryField()[iPatch]).patchNeighbourField();
        }
        else
        {
            psin[iPatch] = sf.boundaryField()[iPatch] + 
                sf.boundaryField()[iPatch].snGrad()
                *
                bmvON_[iPatch]*0.5;
        }
        
        scalarField psio (sf.boundaryField()[iPatch].patchInternalField());
        
        //quad faces
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqx_,0,0,SCA_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqy_,0,1,SCA_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqz_,0,2,SCA_CMPT,VEC_CMPT) //Z
        
        //tri faces
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,batx_,0,0,SCA_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,baty_,0,1,SCA_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,batz_,0,2,SCA_CMPT,VEC_CMPT) //Z
        
        //for other faces - apply surface normal derivative
        {
            label facei = -1;
            forAll(bof_[iPatch], i)
            {
                facei = bof_[iPatch][i];
                gradf.boundaryFieldRef()[iPatch][facei] = 
                    dfdn.boundaryField()[iPatch][facei];
            }
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcGradfIF
(
    const volVectorField& vf,
    const pointVectorField& pf,
    const faceList& faces,
    surfaceTensorField& gradf,
    const surfaceTensorField& dfdn
)
{
    //quad faces
    dfdxif(vf,pf,gradf,qf_,vq_,aqx_,0,0,VEC_CMPT,VEC_CMPT) //X
    dfdxif(vf,pf,gradf,qf_,vq_,aqx_,1,1,VEC_CMPT,VEC_CMPT) //Y
    dfdxif(vf,pf,gradf,qf_,vq_,aqx_,2,2,VEC_CMPT,VEC_CMPT) //Z
    
    dfdxif(vf,pf,gradf,qf_,vq_,aqy_,0,3,VEC_CMPT,VEC_CMPT) //X
    dfdxif(vf,pf,gradf,qf_,vq_,aqy_,1,4,VEC_CMPT,VEC_CMPT) //Y
    dfdxif(vf,pf,gradf,qf_,vq_,aqy_,2,5,VEC_CMPT,VEC_CMPT) //Z
    
    dfdxif(vf,pf,gradf,qf_,vq_,aqz_,0,6,VEC_CMPT,VEC_CMPT) //X
    dfdxif(vf,pf,gradf,qf_,vq_,aqz_,1,7,VEC_CMPT,VEC_CMPT) //Y
    dfdxif(vf,pf,gradf,qf_,vq_,aqz_,2,8,VEC_CMPT,VEC_CMPT) //Z

    //triangular faces
    dfdxif(vf,pf,gradf,tf_,vt_,atx_,0,0,VEC_CMPT,VEC_CMPT) // X
    dfdxif(vf,pf,gradf,tf_,vt_,aty_,1,1,VEC_CMPT,VEC_CMPT) // Y
    dfdxif(vf,pf,gradf,tf_,vt_,atz_,2,2,VEC_CMPT,VEC_CMPT) // Z
    
    dfdxif(vf,pf,gradf,tf_,vt_,atx_,0,3,VEC_CMPT,VEC_CMPT) // X
    dfdxif(vf,pf,gradf,tf_,vt_,aty_,1,4,VEC_CMPT,VEC_CMPT) // Y
    dfdxif(vf,pf,gradf,tf_,vt_,atz_,2,5,VEC_CMPT,VEC_CMPT) // Z
    
    dfdxif(vf,pf,gradf,tf_,vt_,atx_,0,6,VEC_CMPT,VEC_CMPT) // X
    dfdxif(vf,pf,gradf,tf_,vt_,aty_,1,7,VEC_CMPT,VEC_CMPT) // Y
    dfdxif(vf,pf,gradf,tf_,vt_,atz_,2,8,VEC_CMPT,VEC_CMPT) // Z

    //other faces
    {
        label facei = -1;
        forAll(of_, i)
        {
            facei = of_[i];
            gradf.primitiveFieldRef()[facei] = 
                dfdn.primitiveField()[facei];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcGradfBF
(
    const volVectorField& sf,
    const pointVectorField& pf,
    const faceList& faces,
    surfaceTensorField& gradf,
    const surfaceTensorField& dfdn
)
{
    List<List<vector> > psin (sf.boundaryField().size());
    forAll(psin, iPatch)
    {
        if (processorPatch_[iPatch])
        {
            psin[iPatch] = refCast<const processorFvPatchField<vector> >
                (sf.boundaryField()[iPatch]).patchNeighbourField();
        }
        else
        {
            psin[iPatch] = sf.boundaryField()[iPatch] + 
                sf.boundaryField()[iPatch].snGrad()
                *
                bmvON_[iPatch]*0.5;
        }
        
        vectorField psio (sf.boundaryField()[iPatch].patchInternalField());
        
        //quad faces
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqx_,0,0,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqx_,1,1,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqx_,2,2,VEC_CMPT,VEC_CMPT) //Z

        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqy_,0,3,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqy_,1,4,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqy_,2,5,VEC_CMPT,VEC_CMPT) //Z

        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqz_,0,6,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqz_,1,7,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,bqf_,bvq_,baqz_,2,8,VEC_CMPT,VEC_CMPT) //Z

        //tri faces
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,batx_,0,0,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,batx_,1,1,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,batx_,2,2,VEC_CMPT,VEC_CMPT) //Z

        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,baty_,0,3,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,baty_,1,4,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,baty_,2,5,VEC_CMPT,VEC_CMPT) //Z

        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,batz_,0,6,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,batz_,1,7,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,btf_,bvt_,batz_,2,8,VEC_CMPT,VEC_CMPT) //Z

        //for other faces - apply surface normal derivative
        {
            label facei = -1;
            forAll(bof_[iPatch], i)
            {
                facei = bof_[iPatch][i];
                gradf.boundaryFieldRef()[iPatch][facei] = 
                    dfdn.boundaryField()[iPatch][facei];
            }
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::faceGrad(const volScalarField& sf, surfaceVectorField& gradf)
{
    pointScalarField pF
    (
        volPoint_.interpolate
        (
            sf
        )
    );
    //
    const surfaceVectorField& nf = nfRef_();
    surfaceVectorField dfdn
    (
        nf * fvc::snGrad(sf)
    );
    //
    const faceList& faces = sf.mesh().faces();
    /*
     *
     * Calculate grad at internal faces
     *
     */
    calcGradfIF(sf, pF, faces, gradf, dfdn);
    /*
     *
     * Calculate grad at boundary faces
     *
     */
    calcGradfBF(sf, pF, faces, gradf, dfdn);
};

void Foam::fvsc::GaussVolPointBase3D::faceGrad(const volVectorField& vf, surfaceTensorField& gradf)
{
    pointVectorField pF
    (
        volPoint_.interpolate
        (
            vf
        )
    );
    //
    const surfaceVectorField& nf = nfRef_();
    surfaceTensorField dfdn
    (
        nf * fvc::snGrad(vf)
    );
    //
    const faceList& faces = vf.mesh().faces();
    /*
     *
     * Calculate grad at internal faces
     *
     */
    calcGradfIF(vf, pF, faces, gradf, dfdn);
    /*
     *
     * Calculate grad at boundary faces
     *
     */
    calcGradfBF(vf, pF, faces, gradf, dfdn);
};

void Foam::fvsc::GaussVolPointBase3D::faceDiv(const volVectorField& vf, surfaceScalarField& divf)
{
    pointVectorField pF
    (
        volPoint_.interpolate
        (
            vf
        )
    );
    
    surfaceScalarField divfn (nfRef_() & fvc::snGrad(vf));
    const faceList& faces = vf.mesh().faces();
    /*
     *
     * Calculate grad at internal faces
     *
     */
    calcDivfIF(vf, pF, faces, divf, divfn);
    /*
     *
     * Calculate grad at boundary faces
     *
     */
    calcDivfBF(vf, pF, faces, divf, divfn);
};

void Foam::fvsc::GaussVolPointBase3D::faceDiv(const volTensorField& tf, surfaceVectorField& divf)
{
    pointTensorField pF
    (
        volPoint_.interpolate
        (
            tf
        )
    );
    
    surfaceVectorField divfn (nfRef_() & fvc::snGrad(tf));
    const faceList& faces = tf.mesh().faces();
    /*
     *
     * Calculate grad at internal faces
     *
     */
    calcDivfIF(tf, pF, faces, divf, divfn);
    /*
     *
     * Calculate grad at boundary faces
     *
     */
    calcDivfBF(tf, pF, faces, divf, divfn);
}

//
//END-OF-FILE
//


