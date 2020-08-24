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
    Foam::fvsc::GaussVolPoint::GaussVolPoint2D
Description
    This is a method for approximating derivatives of tangents to a face (2D case). 
    They are further used in the calculation of the QGD terms.
\*---------------------------------------------------------------------------*/
#include "GaussVolPointBase2D.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "coupledFvPatch.H"
#include "processorFvPatch.H"
#include "processorFvPatchFields.H"
#include "emptyFvPatch.H"
#include "wedgeFvPatch.H"
#include "volPointInterpolation.H"

Foam::fvsc::GaussVolPointBase2D::GaussVolPointBase2D(const fvMesh& mesh)
:
    c1_(mesh.neighbour().size(), 0.0),
    c2_(mesh.neighbour().size(), 0.0),
    c3_(mesh.neighbour().size(), 0.0),
    c4_(mesh.neighbour().size(), 0.0),
    mv42_(mesh.neighbour().size(), 1.0),
    mv13_(mesh.neighbour().size(), 1.0),
    ip3_(mesh.neighbour().size(), -1),
    ip1_(mesh.neighbour().size(), -1),
    ic4_(mesh.neighbour().size(), -1),
    ic2_(mesh.neighbour().size(), -1),
    e1_(1, 0, 0),
    e2_(0, 1, 0),
    ie1_(-1),
    ie2_(-1),
    ie3_(-1),
    ordinaryPatches_(0),
    ip3e_(mesh.boundary().size()),
    ip1e_(mesh.boundary().size()),
    ic4e_(mesh.boundary().size()),
    c1e_(mesh.boundary().size()),
    c2e_(mesh.boundary().size()),
    c3e_(mesh.boundary().size()),
    c4e_(mesh.boundary().size()),
    mv42e_(mesh.boundary().size()),
    mv13e_(mesh.boundary().size())
{
    if (mesh.nGeometricD() == 2)
    {
        List<scalar> cosa1(mesh.neighbour().size(), 0.0);
        List<scalar> cosa2(mesh.neighbour().size(), 0.0);
        List<scalar> sina1(mesh.neighbour().size(), 0.0);
        List<scalar> sina2(mesh.neighbour().size(), 0.0);
        List<scalar> den(mesh.neighbour().size(), 1.0);
        List<vector> v42(mesh.neighbour().size(), vector::one);
        List<vector> v13(mesh.neighbour().size(), vector::one);
        
        List<List<scalar> > cosa1e(mesh.boundary().size());
        List<List<scalar> > cosa2e(mesh.boundary().size());
        List<List<scalar> > sina1e(mesh.boundary().size());
        List<List<scalar> > sina2e(mesh.boundary().size());
        List<List<scalar> > dene(mesh.boundary().size());
        
        label ip1=-1, ic2=-1, ip3=-1, ic4=-1;
        
        forAll(mesh.geometricD(), iDir)
        {
            if (mesh.geometricD()[iDir] < 1)
            {
                ie3_ = iDir;
            }
        }
        if (ie3_ == 0)
        {
            e1_ = vector(0, 1, 0);
            e2_ = vector(0, 0, 1);
            ie1_ = 1;
            ie2_ = 2;
            ie3_ = 0;
        }
        if (ie3_ == 1)
        {
            e1_ = vector(1, 0, 0);
            e2_ = vector(0, 0, 1);
            ie1_ = 0;
            ie2_ = 2;
            ie3_ = 1;
        }
        if (ie3_ == 2)
        {
            e1_ = vector(1, 0, 0);
            e2_ = vector(0, 1, 0);
            ie1_ = 0;
            ie2_ = 1;
            ie3_ = 2;
        }
        
        forAll(mesh.neighbour(), iFace) //---> for(iFace=0; iFace < mesh_.neighbour().size(); iFace++)
        {
            ic2 = mesh.neighbour()[iFace]; //cell center for point #2
            ic4 = mesh.owner()[iFace];     //cell center for point #4
            
            const face& f = mesh.faces()[iFace];
            #warning "Raise error if number of points on face is not equal to 4"
            forAll(f, ip)
            {
                if (mesh.points()[f[ip]][ie3_] >= mesh.C()[ic2][ie3_])
                {
                    ip1 = f[ip];
                    break;
                }
            }
            forAll(f, ip)
            {
                if (mesh.points()[f[ip]][ie3_] >= mesh.C()[ic2][ie3_])
                {
                    if (ip1 != f[ip])
                    {
                        ip3 = f[ip];
                        break;
                    }
                }
            }
            
            ic2_[iFace]  = ic2;
            ic4_[iFace]  = ic4;
            ip3_[iFace]  = ip3;
            ip1_[iFace]  = ip1;
                         
            v42[iFace]   = mesh.C()[ic2] - mesh.C()[ic4];
            v13[iFace]   = mesh.points()[ip3] - mesh.points()[ip1];
            mv42_[iFace] = mag(v42[iFace]);
            mv13_[iFace] = mag(v13[iFace]);
                         
            cosa1[iFace] = (v42[iFace]/mv42_[iFace]) & e1_;
            cosa2[iFace] = (v13[iFace]/mv13_[iFace]) & e1_;
            sina1[iFace] = (v42[iFace]/mv42_[iFace]) & e2_;
            sina2[iFace] = (v13[iFace]/mv13_[iFace]) & e2_;
                         
            den[iFace]   = sina2[iFace]*cosa1[iFace] - sina1[iFace]*cosa2[iFace];
            c1_[iFace]   = sina2[iFace] / den[iFace];
            c2_[iFace]   = sina1[iFace] / den[iFace];
            c3_[iFace]   = cosa1[iFace] / den[iFace];
            c4_[iFace]   = cosa2[iFace] / den[iFace];
        }
        
        //Info << "Creating weights" << endl;
        forAll(mesh.boundary(), iPatch)
        {
            label patchId = iPatch;
            if (
                !isA<emptyFvPatch>(mesh.boundary()[patchId])
                &&
                !isA<wedgeFvPatch>(mesh.boundary()[patchId])
                )
            {
                if 
                (
                    isA<coupledFvPatch>(mesh.boundary()[patchId])
                    &&
                    !isA<processorFvPatch>(mesh.boundary()[patchId])
                )
                {
                    continue;
                }
                else
                {
                    if (isA<processorFvPatch>(mesh.boundary()[patchId]))
                    {
                        processorPatch_.append(true);
                    }
                    else
                    {
                        processorPatch_.append(false);
                    }
                }
                ordinaryPatches_.append(patchId);
                List<vector> v42(mesh.boundary()[patchId].size(), vector::one);
                List<vector> v13(mesh.boundary()[patchId].size(), vector::one);
                ic4e_[patchId].resize(v42.size());
                mv42e_[patchId].resize(v42.size());
                mv13e_[patchId].resize(v42.size());
                ip3e_[patchId].resize(v42.size());
                ip1e_[patchId].resize(v42.size());
                sina1e[patchId].resize(v42.size());
                sina2e[patchId].resize(v42.size());
                cosa1e[patchId].resize(v42.size());
                cosa2e[patchId].resize(v42.size());
                dene[patchId].resize(v42.size());
                c1e_[patchId].resize(v42.size());
                c2e_[patchId].resize(v42.size());
                c3e_[patchId].resize(v42.size());
                c4e_[patchId].resize(v42.size());
                
                ic4e_[patchId] = mesh.boundary()[patchId].faceCells();
                
                //set v42 vector centers for point 2 if patch is processor
                if (isA<processorFvPatch>(mesh.boundary()[patchId]))
                {
                    const processorPolyPatch & procPolyPatch = 
                        refCast<const processorFvPatch>(mesh.boundary()[patchId]).procPolyPatch();
                    
                    forAll(v42, iFace)
                    {
                        v42[iFace] = procPolyPatch.neighbFaceCellCentres()[iFace] - 
                            mesh.C()[ic4e_[patchId][iFace]];
                    }
                }
                else
                {
                    forAll(v42, iFace)
                    {
                        v42[iFace] = 2.0*(mesh.boundary()[patchId].Cf()[iFace] - 
                            mesh.C()[ic4e_[patchId][iFace]]);
                    }
                }
                
                forAll(mesh.boundary()[patchId], iFace)
                {
                    label ip1=-1, ip3=-1, ic4=ic4e_[patchId][iFace];
                    label globalFaceID = mesh.boundary()[patchId].start() + iFace;
                    
                    const face& f = mesh.faces()[globalFaceID];
                    forAll(f, ip)
                    {
                         if (mesh.points()[f[ip]][ie3_] >= mesh.C()[ic4][ie3_])
                         {
                            ip1 = f[ip];
                            break;
                         }
                    }
                    forAll(f, ip)
                    {
                        if (mesh.points()[f[ip]][ie3_] >= mesh.C()[ic4][ie3_])
                        {
                            if (ip1 != f[ip])
                            {
                                ip3 = f[ip];
                                break;
                            }
                        }
                    }
                    
                    ip3e_[patchId][iFace]  = ip3;
                    ip1e_[patchId][iFace]  = ip1;
                    v13[iFace] = mesh.points()[ip3] - mesh.points()[ip1];
                    
                    mv42e_[patchId][iFace] = mag(v42[iFace]);
                    mv13e_[patchId][iFace] = mag(v13[iFace]);
                    
                    cosa1e[patchId][iFace] = (v42[iFace]/mv42e_[patchId][iFace]) & e1_;
                    cosa2e[patchId][iFace] = (v13[iFace]/mv13e_[patchId][iFace]) & e1_;
                    sina1e[patchId][iFace] = (v42[iFace]/mv42e_[patchId][iFace]) & e2_;
                    sina2e[patchId][iFace] = (v13[iFace]/mv13e_[patchId][iFace]) & e2_;
                    
                    dene[patchId][iFace]   = 
                        sina2e[patchId][iFace]*cosa1e[patchId][iFace]
                        -
                        sina1e[patchId][iFace]*cosa2e[patchId][iFace];
                    
                    c1e_[patchId][iFace]    = sina2e[patchId][iFace] / dene[patchId][iFace];
                    c2e_[patchId][iFace]    = sina1e[patchId][iFace] / dene[patchId][iFace];
                    c3e_[patchId][iFace]    = cosa1e[patchId][iFace] / dene[patchId][iFace];
                    c4e_[patchId][iFace]    = cosa2e[patchId][iFace] / dene[patchId][iFace];
                }
            }
        }
    }
};


Foam::fvsc::GaussVolPointBase2D::~GaussVolPointBase2D()
{
}


void Foam::fvsc::GaussVolPointBase2D::faceGrad(const volScalarField& f, surfaceVectorField& gradf)
{
    scalar dfdn = 0.0, dfdt = 0.0;
    
    if (f.mesh().nGeometricD() == 2)
    {
        pointScalarField pF
        (
            volPointInterpolation::New(f.mesh()).interpolate
            (
                f
            )
        );
        
        forAll(f.mesh().neighbour(), iFace) //---> for(iFace=0; iFace < mesh_.neighbour().size(); iFace++)
        {
            dfdn = (f[ic2_[iFace]] - f[ic4_[iFace]]) / mv42_[iFace];
            //
            
            // df/dl2 (1-3)
            dfdt = (pF[ip3_[iFace]] - pF[ip1_[iFace]]) / mv13_[iFace];
            //
            
            gradf[iFace][ie1_] =  (dfdn*c1_[iFace] - dfdt*c2_[iFace]);
            gradf[iFace][ie2_] =  (dfdt*c3_[iFace] - dfdn*c4_[iFace]);
            //gradf[iFace][ie1_] =  (dfdn*sina2_[iFace] - dfdt*sina1_[iFace])/den_[iFace];
            //gradf[iFace][ie2_] =  (dfdt*cosa1_[iFace] - dfdn*cosa2_[iFace])/den_[iFace];
            gradf[iFace][ie3_] = 0.0;
        }
        
        List<List<scalar> > psi2(f.boundaryField().size());
        
        forAll(ordinaryPatches_, iPatch)
        {
            label patchId = ordinaryPatches_[iPatch];
            if (processorPatch_[iPatch])
            {
                psi2[patchId] = refCast<const processorFvPatchField<scalar> >
                    (f.boundaryField()[patchId]).patchNeighbourField();
            }
            else
            {
                psi2[patchId] = f.boundaryField()[patchId] + 
                    f.boundaryField()[patchId].snGrad()
                    *
                    mv42e_[patchId]*0.5;
            }
            
            forAll(f.boundaryField()[patchId], iFace)
            {
                //dfdn
                dfdn = (psi2[patchId][iFace] - f[ic4e_[patchId][iFace]]) / mv42e_[patchId][iFace];
                
                //dfdt
                dfdt = (pF[ip3e_[patchId][iFace]] - pF[ip1e_[patchId][iFace]]) / mv13e_[patchId][iFace];
                
                gradf.boundaryFieldRef()[patchId][iFace][ie1_] = (dfdn*c1e_[patchId][iFace] - dfdt*c2e_[patchId][iFace]);
                gradf.boundaryFieldRef()[patchId][iFace][ie2_] = (dfdt*c3e_[patchId][iFace] - dfdn*c4e_[patchId][iFace]);
                gradf.boundaryFieldRef()[patchId][iFace][ie3_] = 0.0;
            }
        }
    }
    else
    {
        #warning "Raise error if called not for 2D case"
    }
};

void Foam::fvsc::GaussVolPointBase2D::faceDiv(const volVectorField& f, surfaceScalarField& divf)
{
    scalar df1dn = 0.0, df1dt = 0.0,
        df2dn = 0.0, df2dt = 0.0;
    
    if (f.mesh().nGeometricD() == 2)
    {
        pointVectorField pF
        (
            volPointInterpolation::New(f.mesh()).interpolate
            (
                f
            )
        );
        
        forAll(f.mesh().neighbour(), iFace) //---> for(iFace=0; iFace < mesh_.neighbour().size(); iFace++)
        {
            df1dn = (f[ic2_[iFace]][ie1_] - f[ic4_[iFace]][ie1_]) / mv42_[iFace];
            df2dn = (f[ic2_[iFace]][ie2_] - f[ic4_[iFace]][ie2_]) / mv42_[iFace];
            //
            
            // df/dl2 (1-3)
            df1dt = (pF[ip3_[iFace]][ie1_] - pF[ip1_[iFace]][ie1_]) / mv13_[iFace];
            df2dt = (pF[ip3_[iFace]][ie2_] - pF[ip1_[iFace]][ie2_]) / mv13_[iFace];
            //
            
            divf[iFace] = (df1dn*c1_[iFace] - df1dt*c2_[iFace]) + 
                (df2dt*c3_[iFace] - df2dn*c4_[iFace]);
        }
        
        List<List<vector> > psi2(f.boundaryField().size());
        
        forAll(ordinaryPatches_, iPatch)
        {
            label patchId = ordinaryPatches_[iPatch];
            if (processorPatch_[iPatch])
            {
                psi2[patchId] = refCast<const processorFvPatchField<vector> >
                    (f.boundaryField()[patchId]).patchNeighbourField();
            }
            else
            {
                psi2[patchId] = f.boundaryField()[patchId] + 
                    f.boundaryField()[patchId].snGrad()
                    *
                    mv42e_[patchId]*0.5;
            }
            
            forAll(f.boundaryField()[patchId], iFace)
            {
                //dfdn
                df1dn = (psi2[patchId][iFace][ie1_] - f[ic4e_[patchId][iFace]][ie1_]) / mv42e_[patchId][iFace];
                df2dn = (psi2[patchId][iFace][ie2_] - f[ic4e_[patchId][iFace]][ie2_]) / mv42e_[patchId][iFace];
                
                //dfdt
                df1dt = (pF[ip3e_[patchId][iFace]][ie1_] - pF[ip1e_[patchId][iFace]][ie1_]) / mv13e_[patchId][iFace];
                df2dt = (pF[ip3e_[patchId][iFace]][ie2_] - pF[ip1e_[patchId][iFace]][ie2_]) / mv13e_[patchId][iFace];
                
                divf.boundaryFieldRef()[patchId][iFace] = 
                    (df1dn*c1e_[patchId][iFace] - df1dt*c2e_[patchId][iFace])
                    +
                    (df2dt*c3e_[patchId][iFace] - df2dn*c4e_[patchId][iFace]);
            }
        }
    }
    else
    {
        #warning "Raise error if called not for 2D case"
    }
};

void Foam::fvsc::GaussVolPointBase2D::faceDiv(const volTensorField& f, surfaceVectorField& divf)
{
    scalar df11dn = 0.0, df11dt = 0.0,
           df21dn = 0.0, df21dt = 0.0,
           df22dn = 0.0, df22dt = 0.0,
           df12dn = 0.0, df12dt = 0.0;
    
    label i11 = ie1_*3 + ie1_;
    label i21 = ie2_*3 + ie1_;
    label i22 = ie2_*3 + ie2_;
    label i12 = ie1_*3 + ie2_;
    
    if (f.mesh().nGeometricD() == 2)
    {
        pointTensorField pF
        (
            volPointInterpolation::New(f.mesh()).interpolate
            (
                f
            )
        );
        
        forAll(f.mesh().neighbour(), iFace) //---> for(iFace=0; iFace < mesh_.neighbour().size(); iFace++)
        {
            df11dn = (f[ic2_[iFace]][i11] - f[ic4_[iFace]][i11]) / mv42_[iFace];
            df21dn = (f[ic2_[iFace]][i21] - f[ic4_[iFace]][i21]) / mv42_[iFace];
            
            // df/dl2 (1-3)
            df11dt = (pF[ip3_[iFace]][i11] - pF[ip1_[iFace]][i11]) / mv13_[iFace];
            df21dt = (pF[ip3_[iFace]][i21] - pF[ip1_[iFace]][i21]) / mv13_[iFace];
            
            divf[iFace][ie1_] = 
                (df11dn*c1_[iFace] - df11dt*c2_[iFace]) + //d/dx
                (df21dt*c3_[iFace] - df21dn*c4_[iFace]);  //d/dy
            
            df22dn = (f[ic2_[iFace]][i22] - f[ic4_[iFace]][i22]) / mv42_[iFace];
            df12dn = (f[ic2_[iFace]][i12] - f[ic4_[iFace]][i12]) / mv42_[iFace];
            
            // df/dl2 (1-3)
            df22dt = (pF[ip3_[iFace]][i22] - pF[ip1_[iFace]][i22]) / mv13_[iFace];
            df12dt = (pF[ip3_[iFace]][i12] - pF[ip1_[iFace]][i12]) / mv13_[iFace];

            divf[iFace][ie2_] = 
                (df12dn*c1_[iFace] - df12dt*c2_[iFace]) + //d/dx
                (df22dt*c3_[iFace] - df22dn*c4_[iFace]);  //d/dy
        }
        
        List<List<tensor> > psi2(f.boundaryField().size());
        
        forAll(ordinaryPatches_, iPatch)
        {
            label patchId = ordinaryPatches_[iPatch];
            if (processorPatch_[iPatch])
            {
                psi2[patchId] = refCast<const processorFvPatchField<tensor> >
                    (f.boundaryField()[patchId]).patchNeighbourField();
            }
            else
            {
                psi2[patchId] = f.boundaryField()[patchId] + 
                    f.boundaryField()[patchId].snGrad()
                    *
                    mv42e_[patchId]*0.5;
            }
            
            forAll(f.boundaryField()[patchId], iFace)
            {
                //dfdn
                df11dn = (psi2[patchId][iFace][i11] - f[ic4e_[patchId][iFace]][i11]) / mv42e_[patchId][iFace];
                df21dn = (psi2[patchId][iFace][i21] - f[ic4e_[patchId][iFace]][i21]) / mv42e_[patchId][iFace];
                
                //dfdt
                df11dt = (pF[ip3e_[patchId][iFace]][i11] - pF[ip1e_[patchId][iFace]][i11]) / mv13e_[patchId][iFace];
                df21dt = (pF[ip3e_[patchId][iFace]][i21] - pF[ip1e_[patchId][iFace]][i21]) / mv13e_[patchId][iFace];
                
                divf.boundaryFieldRef()[patchId][iFace][ie1_] = 
                    (df11dn*c1e_[patchId][iFace] - df11dt*c2e_[patchId][iFace])  //d/dx
                    +
                    (df21dt*c3e_[patchId][iFace] - df21dn*c4e_[patchId][iFace]); //d/dy
                
                //dfdn
                df12dn = (psi2[patchId][iFace][i12] - f[ic4e_[patchId][iFace]][i12]) / mv42e_[patchId][iFace];
                df22dn = (psi2[patchId][iFace][i22] - f[ic4e_[patchId][iFace]][i22]) / mv42e_[patchId][iFace];
                
                //dfdt
                df12dt = (pF[ip3e_[patchId][iFace]][i12] - pF[ip1e_[patchId][iFace]][i12]) / mv13e_[patchId][iFace];
                df22dt = (pF[ip3e_[patchId][iFace]][i22] - pF[ip1e_[patchId][iFace]][i22]) / mv13e_[patchId][iFace];
                
                divf.boundaryFieldRef()[patchId][iFace][ie2_] = 
                    (df12dn*c1e_[patchId][iFace] - df12dt*c2e_[patchId][iFace])  //d/dx
                    +
                    (df22dt*c3e_[patchId][iFace] - df22dn*c4e_[patchId][iFace]); //d/dy
            }
        }
    }
    else
    {
        #warning "Raise error if called not for 2D case"
    }
}

//
//END-OF-FILE
//


