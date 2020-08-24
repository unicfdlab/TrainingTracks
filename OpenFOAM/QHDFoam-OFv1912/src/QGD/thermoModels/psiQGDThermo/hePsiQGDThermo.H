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
    Foam::hePsiQGDThermo

Group
    grpPsiQGDThermo

Description
    Energy for a perfect gas mixture with QGD equations

SourceFiles
    hePsiQGDThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hePsiQGDThermo_H
#define hePsiQGDThermo_H

#include "psiThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class hePsiQGDThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class hePsiQGDThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        hePsiQGDThermo(const hePsiQGDThermo<BasicPsiThermo, MixtureType>&);

public:

    //- Runtime type information
    TypeName("hePsiQGDThermo");


    // Constructors

        //- Construct from mesh and phase name
        hePsiQGDThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and phase name
        hePsiQGDThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );


    //- Destructor
    virtual ~hePsiQGDThermo();


    // Member functions

        //- Update properties
        virtual void correct();

    // Accesss functions

    virtual tmp<volScalarField> gamma() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "hePsiQGDThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
