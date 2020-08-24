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
    Foam::rhoQGDThermo

Group
    grpRhoQGDThermo

Description
    Basic thermodynamic properties based on density

SourceFiles
    rhoQGDThermo.C

\*---------------------------------------------------------------------------*/

#ifndef rhoQGDThermo_H
#define rhoQGDThermo_H

#include "rhoThermo.H"
#include "QGDThermo.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class rhoQGDThermo Declaration
\*---------------------------------------------------------------------------*/

class rhoQGDThermo
:
    public rhoThermo, public QGDThermo
{

protected:

    // Protected Member Functions

        //- Construct as copy (not implemented)
        rhoQGDThermo(const rhoQGDThermo&);

        //-
        volScalarField c_;

        //-
        volScalarField gamma_;


public:

    //- Runtime type information
    TypeName("rhoQGDThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        rhoQGDThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        rhoQGDThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and phase name
        rhoQGDThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );

    //- Selector
    static autoPtr<rhoQGDThermo> New
    (
        const fvMesh& mesh,
        const word& phaseName=word::null
    );

    //- Selector
    static autoPtr<rhoQGDThermo> New
    (
        const fvMesh& mesh,
        const word& phaseName,
        const word& dictName
    );


    //- Destructor
    virtual ~rhoQGDThermo();


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
