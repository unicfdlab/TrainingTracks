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
    Foam::powerLawTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    powerLawTransportI.H
    powerLawTransport.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawTransport_H
#define powerLawTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class powerLawTransport;

template<class Thermo>
inline powerLawTransport<Thermo> operator+
(
    const powerLawTransport<Thermo>&,
    const powerLawTransport<Thermo>&
);

template<class Thermo>
inline powerLawTransport<Thermo> operator-
(
    const powerLawTransport<Thermo>&,
    const powerLawTransport<Thermo>&
);

template<class Thermo>
inline powerLawTransport<Thermo> operator*
(
    const scalar,
    const powerLawTransport<Thermo>&
);

template<class Thermo>
inline powerLawTransport<Thermo> operator==
(
    const powerLawTransport<Thermo>&,
    const powerLawTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const powerLawTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class powerLawTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class powerLawTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa.s]
        scalar mu0_;

        //- [K]
        scalar T0_;
        
        //- power coeff.
        scalar k_;

        //- Reciprocal Prandtl Number []
        scalar rPr_;


    // Private Member Functions

        //- Construct from components
        inline powerLawTransport
        (
            const Thermo& t,
            const scalar mu0,
            const scalar T0,
            const scalar kk,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline powerLawTransport(const word&, const powerLawTransport&);

        //- Construct from Istream
        powerLawTransport(Istream&);

        //- Construct from dictionary
        powerLawTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<powerLawTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<powerLawTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<powerLawTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "powerLaw<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline powerLawTransport& operator=(const powerLawTransport&);

        inline void operator+=(const powerLawTransport&);

        inline void operator-=(const powerLawTransport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend powerLawTransport operator+ <Thermo>
        (
            const powerLawTransport&,
            const powerLawTransport&
        );

        friend powerLawTransport operator- <Thermo>
        (
            const powerLawTransport&,
            const powerLawTransport&
        );

        friend powerLawTransport operator* <Thermo>
        (
            const scalar,
            const powerLawTransport&
        );

        friend powerLawTransport operator== <Thermo>
        (
            const powerLawTransport&,
            const powerLawTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const powerLawTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "powerLawTransportI.H"

#ifdef NoRepository
#   include "powerLawTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
