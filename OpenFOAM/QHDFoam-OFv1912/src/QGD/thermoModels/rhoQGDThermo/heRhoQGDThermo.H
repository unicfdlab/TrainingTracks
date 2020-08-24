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
    Foam::heRhoThermo

Group
    grpRhoQGDThermo

Description
    Energy for a mixture based on density

SourceFiles
    heRhoQGDThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heRhoQGDThermo_H
#define heRhoQGDThermo_H

#include "rhoThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heRhoQGDThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class heRhoQGDThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate();

        //- Construct as copy (not implemented)
        heRhoQGDThermo(const heRhoQGDThermo<BasicPsiThermo, MixtureType>&);

public:

    //- Runtime type information
    TypeName("heRhoQGDThermo");


    // Constructors

        //- Construct from mesh and phase name
        heRhoQGDThermo
        (
            const fvMesh&,
            const word& phaseName
        );


        //- Construct from mesh and phase name and dictionary name
        heRhoQGDThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );

    //- Destructor
    virtual ~heRhoQGDThermo();


    // Member functions

        //- Update properties
        virtual void correct();
        
    // Accesss functions

    //virtual tmp<volScalarField> gamma() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "heRhoQGDThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
