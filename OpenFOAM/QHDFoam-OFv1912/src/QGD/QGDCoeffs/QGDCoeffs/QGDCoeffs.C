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

#include "QGDCoeffs.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledFvsPatchFields.H"
#include "QGDThermo.H"
#include "linear.H"
#include "emptyFvPatch.H"
#include "wedgeFvPatch.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(QGDCoeffs, 0);
    defineRunTimeSelectionTable(QGDCoeffs, dictionary);
}
}

namespace Foam
{
namespace qgd
{

autoPtr<QGDCoeffs> QGDCoeffs::New
(
    const word& qgdCoeffsType,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    Info<< "Selecting QGD coeffs evaluation approach type " << qgdCoeffsType << endl;
    
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(qgdCoeffsType);
    
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "QGDCoeffs::New(const word&, const fvMesh&)"
        )   << "Unknown QGD coeffs evaluation approach type " << qgdCoeffsType << nl << nl
        << "Valid model types are:" << nl
        << dictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }
    
    if (dict.found(qgdCoeffsType + "Dict"))
    {
        return autoPtr<QGDCoeffs>
        (
            cstrIter()
            (
                IOobject
                (
                    qgdCoeffsType,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dict.subDict(qgdCoeffsType + "Dict")
            )
         );
    }
    
    return autoPtr<QGDCoeffs>
    (
        cstrIter()
        (
            IOobject
            (
                qgdCoeffsType,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dict
        )
    );
}

tmp<volScalarField> QGDCoeffs::readOrCreateAlphaQGD(const fvMesh& mesh)
{
    IOobject aQGDHeader
    (
        "alphaQGD",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );
    
    if (aQGDHeader.typeHeaderOk<volScalarField>())
    {
        aQGDHeader.writeOpt() = IOobject::AUTO_WRITE;
        
        return
        tmp<volScalarField>
        (
            new volScalarField
            (
                aQGDHeader,
                mesh
            )
        );
    }
    
    tmp<volScalarField> newAlphaQGD
    (
        new volScalarField
        (
            aQGDHeader,
            mesh,
            dimensionSet(0, 0, 0, 0, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    
    newAlphaQGD.ref().primitiveFieldRef() = 0.5;
    newAlphaQGD.ref().correctBoundaryConditions();
    
    return newAlphaQGD;
}

//
QGDCoeffs::QGDCoeffs(const IOobject& io, const fvMesh& mesh, const dictionary& dict)
:
    regIOobject(io, false),
    refCount(),
    mesh_(mesh),
    runTime_(mesh_.time()),
    muQGD_
    (
        IOobject
        (
            "muQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),
    alphauQGD_
    (
        IOobject
        (
            "alphauQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),
    hQGDf_
    (
        "hQGDf",
        1.0 / mag(mesh.surfaceInterpolation::deltaCoeffs())
    ),
    hQGD_
    (
        IOobject
        (
            "hQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        hQGDf_.dimensions()
    ),
    taQGD_
    (
        readOrCreateAlphaQGD(mesh)
    ),
    aQGD_(taQGD_.ref()),
    tauQGD_
    (
        IOobject
        (
            "tauQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 1, 0, 0)
    ),
    tauQGDf_
    (
        "tauQGDf",
        linearInterpolate(tauQGD_)
    ),
    PrQGD_
    (
        IOobject
        (
            "PrQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    ),
    ScQGD_
    (
        IOobject
        (
            "ScQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    )
{
    this->updateQGDLength(mesh);
}

QGDCoeffs::~QGDCoeffs()
{
}

void Foam::qgd::QGDCoeffs::correct(const QGDThermo& qgdThermo)
{
    forAll(tauQGD_, celli)
    {
        tauQGD_.primitiveFieldRef()[celli] = 0.0;
        muQGD_.primitiveFieldRef()[celli] = 0.0;
        alphauQGD_.primitiveFieldRef()[celli] = 0.0;
        ScQGD_.primitiveFieldRef()[celli] = 1.0;
        PrQGD_.primitiveFieldRef()[celli] = 1.0;
    }
    forAll(tauQGD_.boundaryField(), patchi)
    {
        forAll(tauQGD_.boundaryField()[patchi], facei)
        {
            tauQGD_.boundaryFieldRef()[patchi][facei] =
                0.0;
            muQGD_.boundaryFieldRef()[patchi][facei] =
                0.0;
            alphauQGD_.boundaryFieldRef()[patchi][facei] =
                0.0;
            PrQGD_.boundaryFieldRef()[patchi][facei] =
                1.0;
            ScQGD_.boundaryFieldRef()[patchi][facei] =
                1.0;
        }
    }
}

void Foam::qgd::QGDCoeffs::updateQGDLength(const fvMesh& mesh)
{
    {
        scalar hown = 0.0;
        scalar hnei = 0.0;
        forAll(hQGDf_.primitiveField(), iFace)
        {
            hown = mag(mesh.C()[mesh.owner()[iFace]] - mesh.Cf()[iFace]);
            hnei = mag(mesh.C()[mesh.neighbour()[iFace]] - mesh.Cf()[iFace]);
            hQGDf_.primitiveFieldRef()[iFace] = 2.0*min(hown, hnei);
        }

        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            if (!fvp.coupled())
            {
                hQGDf_.boundaryFieldRef()[patchi] *= 2.0;
            }
        }
    }

    scalar hint = 0.0;
    scalar surf = 0.0;
    label  fid  = 0;
    forAll(hQGD_, celli)
    {
        const cell& c = mesh.cells()[celli];
        hint = 0.0;
        surf = 0.0;
        label pid = -1;
        label pfid = -1;
        forAll(c, facei)
        {
            fid  = c[facei];

            if (mesh.isInternalFace(fid))
            {
                hint += hQGDf_[fid] * mesh.magSf()[fid];
                surf += mesh.magSf()[fid];
            }
            else
            {
                pid = mesh.boundaryMesh().whichPatch(fid);
                if (pid >= 0)
                {
                    if 
                    (
                        !isA<emptyFvPatch>(mesh.boundary()[pid])
                        &&
                        !isA<wedgeFvPatch>(mesh.boundary()[pid])
                    )
                    {
                        pfid = mesh.boundaryMesh()[pid].whichFace(fid);

                        hint += hQGDf_.boundaryField()[pid][pfid] *
                            mesh.magSf().boundaryField()[pid][pfid];
                        surf += mesh.magSf().boundaryField()[pid][pfid];
                    }
                }
            }
        }

        hQGD_.primitiveFieldRef()[celli] = hint / surf;
    }

    forAll(mesh.boundary(), patchi)
    {
        if 
        (
        !isA<emptyFvPatch>(mesh.boundary()[patchi])
        //&&
        //!isA<wedgeFvPatch>(mesh.boundary()[patchi])
        )
        {
            hQGD_.boundaryFieldRef()[patchi] = hQGDf_.boundaryField()[patchi] * 1.0;
        }
    }
}

const Foam::surfaceScalarField& Foam::qgd::QGDCoeffs::hQGDf() const
{
    return hQGDf_;
}

const Foam::volScalarField& Foam::qgd::QGDCoeffs::hQGD() const
{
    return hQGD_;
}

const Foam::volScalarField& Foam::qgd::QGDCoeffs::alphauQGD() const
{
    return alphauQGD_;
}

const Foam::volScalarField& Foam::qgd::QGDCoeffs::muQGD() const
{
    return muQGD_;
}

const Foam::volScalarField& Foam::qgd::QGDCoeffs::tauQGD() const
{
    return tauQGD_;
}

const Foam::surfaceScalarField& Foam::qgd::QGDCoeffs::tauQGDf() const
{
    return tauQGDf_;
}

}; //namespace qgd

}; //namespace Foam


//END-OF-FILE
