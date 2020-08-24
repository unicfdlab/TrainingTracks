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

Class
    Foam::qgd::QGDCoeffs

Group
    grpQGDCoeffs

Description
    Base class for all classes describing QGD model coefficients.
    Provides interfaces for accessing QGD/QHD models coefficients \mu^{QGD},
    \tau^{QGD}, \alpha^{QGD}, \alpha_u^{QGD}, h^{QGD}

SourceFiles
    QGDCoeffs.C

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "regIOobject.H"
#include "volFields.H"
#include "surfaceFields.H"

#ifndef QGDCoeffs_H
#define QGDCoeffs_H

namespace Foam
{

//forward declaration of psiQGDThermo
class QGDThermo;

namespace qgd
{

class QGDCoeffs : public regIOobject, public refCount
{

protected:

    //-
    const fvMesh& mesh_;

    //-
    const Time& runTime_;

    //-
    volScalarField muQGD_;

    //-
    volScalarField alphauQGD_;

    //-
    surfaceScalarField hQGDf_;

    //-
    volScalarField hQGD_;

    //-
    tmp<volScalarField> taQGD_;
    
    //-
    volScalarField& aQGD_;

    //-
    volScalarField tauQGD_;

    //-
    surfaceScalarField tauQGDf_;

    //-
    volScalarField PrQGD_;

    //-
    volScalarField ScQGD_;

protected:

    //-
    virtual void updateQGDLength(const fvMesh&);
    
    //-
    tmp<volScalarField> readOrCreateAlphaQGD(const fvMesh& );

public:

    //-
    TypeName("QGDCoeffs");

    //-
    declareRunTimeSelectionTable
    (
        autoPtr,
        QGDCoeffs,
        dictionary,
        (
            const IOobject& io,
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (io, mesh, dict)
    );

    //- Construct from components
    QGDCoeffs
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //- Return a reference to the selected fvscStencil model
    static autoPtr<QGDCoeffs> New
    (
        const word& qgdCoeffsType,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    virtual ~QGDCoeffs();

    //-
    virtual bool writeData(Ostream&) const
    {
        return true;
    }

    //-
    virtual void correct(const QGDThermo&);

    //-
    const volScalarField& hQGD() const;

    //-
    const volScalarField& tauQGD() const;

    //-
    const surfaceScalarField& hQGDf() const;

    //-
    const surfaceScalarField& tauQGDf() const;

    //-
    const volScalarField& muQGD() const;

    //-
    const volScalarField& alphauQGD() const;

};

}

}

#endif

//
//END-OF-FILE
//
