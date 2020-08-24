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
    This file is part of QGDsolver, based on OpenFOAM library.

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
    Foam::psiQGDThermo

Group
    grpPsiQGDThermo

Description
    Class describing thermophysical properties of perfect gas with motion
    governed by Quasi- Gas dynamics equations.

SourceFiles
    psiQGDThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiQGDThermo_H
#define psiQGDThermo_H

#include "QGDThermo.H"
#include "psiThermo.H"
#include "runTimeSelectionTables.H"
#include "QGDCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class psiQGDThermo Declaration
\*---------------------------------------------------------------------------*/

class psiQGDThermo
:
    public psiThermo, public QGDThermo
{

protected:

    // Protected Member Functions

        //- Construct as copy (not implemented)
        psiQGDThermo(const psiQGDThermo&);

        //-
        volScalarField c_;

        //-
        volScalarField gamma_;

public:

    //- Runtime type information
    TypeName("psiQGDThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiQGDThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        psiQGDThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and phase name
        psiQGDThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Selector
    static autoPtr<psiQGDThermo> New
    (
        const fvMesh& mesh,
        const word& phaseName=word::null
    );

    //- Selector
    static autoPtr<psiQGDThermo> New
    (
        const fvMesh& mesh,
        const word& phaseName,
        const word& dictName
    );

    //- Destructor
    virtual ~psiQGDThermo();


    // Member functions

    //-
    bool read();

    //-
    const volScalarField& c() const;

    //-
    virtual const volScalarField& p() const;

    //-
    volScalarField& p();

    //-
    virtual tmp<volScalarField> rho() const;

    //-
    virtual tmp<volScalarField> mu() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
