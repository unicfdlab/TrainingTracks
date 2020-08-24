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
\*---------------------------------------------------------------------------*/

#include "varScModel6.H"
#include "psiQGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSmooth.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "primitiveMeshTools.H"
#include "linear.H"
#include "cellSet.H"
#include "emptyFvPatch.H"
#include "wedgeFvPatch.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(varScModel6,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        varScModel6,
        dictionary
    );
}
}

Foam::qgd::
varScModel6::varScModel6
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict),
    smoothCoeff_(0.1),
    rC_(0.5),
    minSc_(0.05),
    maxSc_(1.00),
    badQualitySc_(0.05),
    qgdAspectRatioThreshold_(1.5),
    constSc_(0.05),
    cqSc_(mesh.V().size(), 0.0),
    constScCellSetPtr_(nullptr)
{
    scalar ScQGD = 0.0, PrQGD = 1.0;

    dict.lookup("ScQGD") >> ScQGD;
    dict.lookup("PrQGD") >> PrQGD;

    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;

    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;

    if (dict.found("smoothCoeff"))
    {
        dict.lookup("smoothCoeff") >> smoothCoeff_;
    }

    if (dict.found("rC"))
    {
        dict.lookup("rC") >> rC_;
    }

    if (dict.found("minSc"))
    {
        dict.lookup("minSc") >> minSc_;
    }

    if (dict.found("maxSc"))
    {
        dict.lookup("maxSc") >> maxSc_;
    }

    if (dict.found("badQualitySc"))
    {
        dict.lookup("badQualitySc") >> badQualitySc_;
    }
    
    if (dict.found("maxAspectRatio"))
    {
        dict.lookup("maxAspectRatio") >> qgdAspectRatioThreshold_;
    }

    scalarField openness(mesh.V().size(), 0);
    scalarField aspectRatio(mesh.V().size(), 1);

    primitiveMeshTools::cellClosedness
    (
        mesh_,
        mesh_.geometricD(),
        mesh_.faceAreas(),
        mesh_.cellVolumes(),
        openness,
        aspectRatio
    );

    forAll(aspectRatio, iCell)
    {
        if (aspectRatio[iCell] > qgdAspectRatioThreshold_)
        {
            cqSc_[iCell] = badQualitySc_ * aspectRatio[iCell] / qgdAspectRatioThreshold_;
            //cqSc_[iCell] = badQualitySc_;
        }
    }

    if (dict.found("constScCellSet"))
    {
        word constScCellSet (dict.lookup("constScCellSet"));
        constSc_ = ScQGD;

        constScCellSetPtr_.reset
        (
            new cellSet
            (
                mesh,
                constScCellSet,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
    }

    /*
    {
        volScalarField hRatio(this->hQGD_);
        label fid = -1, pid = -1, pfid = -1;
        scalar hmax = 0.0;
        forAll(hRatio, celli)
        {
            const cell& c = mesh.cells()[celli];
            hmax = 0.0;
            forAll(c, facei)
            {
                fid  = c[facei];
                if (mesh.isInternalFace(fid))
                {
                    if (hQGDf_[fid] > hmax)
                    {
                        hmax = hQGDf_[fid];
                    }
                }
                else
                {
                    pid = mesh.boundaryMesh().whichPatch(fid);
                    pfid= mesh.boundary()[pid].patch().whichFace(fid);

                    if (hQGDf_.boundaryField()[pid][pfid] > hmax)
                    {
                        hmax = hQGDf_.boundaryField()[pid][pfid];
                    }
                }
            }
            hRatio.primitiveFieldRef()[celli] = hmax / hQGD_[celli];
            if (hRatio[celli] > 1.2)
            {
                cqSc_[celli] = 1.0;
            }
        }
    }
    */
    //fvc::smooth(this->hQGD_, smoothCoeff_);
}

Foam::qgd::
varScModel6::~varScModel6()
{
}

void Foam::qgd::
varScModel6::correct(const Foam::QGDThermo& qgdThermo)
{
    const volScalarField& cSound = qgdThermo.c();
    const volScalarField& p      = qgdThermo.p();
    const volScalarField rho(qgdThermo.rho()*1.0);

    this->tauQGDf_= linearInterpolate(this->aQGD_ / cSound) * hQGDf_;
    this->tauQGD_ = this->aQGD_ * this->hQGD_  / cSound;
    
    const surfaceScalarField pf =
        linearInterpolate(p);
    const surfaceScalarField dpf =
        fvc::snGrad(p)/mesh_.deltaCoeffs();
    
    forAll(ScQGD_, icell)
    {
        const cell& facesOfCell = mesh_.cells()[icell];
        ScQGD_[icell] = 0.0;
        scalar sumDpF = 0.0;
        scalar sumpf  = 0.0;
        label faceid = -1;
        label patchid = -1;
        scalar n = 0.0;
        
        forAll(facesOfCell, iface)
        {
            faceid  = facesOfCell[iface];
            if (mesh_.isInternalFace(faceid))
            {
                sumpf  +=  pf.primitiveField()[faceid];
                if (mesh_.owner()[faceid] == icell)
                {
                    sumDpF += dpf.primitiveField()[faceid];
                }
                else
                {
                    sumDpF -= dpf.primitiveField()[faceid];
                }
                n = n + 1;
            }
            else
            {
                patchid  = mesh_.boundaryMesh().whichPatch(faceid);
                if (isA<emptyFvPatch>(mesh_.boundary()[patchid]))
                {
                    continue;
                }
                if (isA<wedgeFvPatch>(mesh_.boundary()[patchid]))
                {
                    continue;
                }
                faceid  -= mesh_.boundaryMesh()[patchid].start();
                
                if (p.boundaryField()[patchid].coupled())
                {
                    sumDpF += 
                        (p.boundaryField()[patchid].patchNeighbourField()()[faceid] - p[icell]);
                }
                else
                {
                    sumDpF += dpf.boundaryField()[patchid][faceid];
                }
                sumpf  +=  pf.boundaryField()[patchid][faceid];
                n = n + 1;
            }
        }
        sumpf /= n;
        ScQGD_[icell] = mag(sumDpF) / sumpf;
    }


//    if (mesh_.time().outputTime())
//    {
//        sumdp.write();
//        avrp.write();
//    }

//    this->ScQGD_ =
//        rC_ *
//        (mag(fvc::grad(rho)) * hQGD_ / rho) +
//        (1.0 - rC_) * ScQGD_;

//    this->ScQGD_ =
//        max(this->ScQGD_, minSc_);
//    this->ScQGD_ =
//        min(this->ScQGD_, maxSc_);
//
//    this->ScQGD_.primitiveFieldRef() =
//        max(this->ScQGD_.primitiveField(), cqSc_);
//
//    if (constScCellSetPtr_.valid())
//    {
//        const cellSet& constScCells = constScCellSetPtr_();
//
//        forAllConstIter(cellSet, constScCells, iter)
//        {
//            ScQGD_[iter.key()] = constSc_;
//        }
//    }
//
//    fvc::smooth(this->ScQGD_, smoothCoeff_);

    Info<< "max/min ScQGD: "
        << max(this->ScQGD_).value() << "/"
        << min(this->ScQGD_).value() << endl;

    if (runTime_.outputTime())
    {
        this->ScQGD_.write();
        this->hQGD_.write();
    }

    forAll(p.primitiveField(), celli)
    {
        muQGD_.primitiveFieldRef()[celli] =
            p.primitiveField()[celli] *
            ScQGD_.primitiveField()[celli] *
            tauQGD_.primitiveField()[celli];

        alphauQGD_.primitiveFieldRef()[celli] = muQGD_.primitiveField()[celli] /
            PrQGD_.primitiveField()[celli];
    }

    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            muQGD_.boundaryFieldRef()[patchi][facei] =
                p.boundaryField()[patchi][facei] *
                ScQGD_.boundaryField()[patchi][facei] *
                tauQGD_.boundaryField()[patchi][facei];

            alphauQGD_.boundaryFieldRef()[patchi][facei] =
                muQGD_.boundaryFieldRef()[patchi][facei] /
                PrQGD_.boundaryField()[patchi][facei];
        }
    }
}

//
//END-OF-FILE
//
