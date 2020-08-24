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

Group
    grpPowerLawQGD

Description:
    Definition of templates

SourceFiles:
    powerLawTransportI.H
    powerLawTransport.C

\*----------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawTransport<Thermo>::powerLawTransport
(
    const Thermo& t,
    const scalar mu0,
    const scalar T0,
    const scalar kk,
    const scalar Pr
)
:
    Thermo(t),
    mu0_(mu0),
    T0_(T0),
    k_(kk),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::powerLawTransport<Thermo>::powerLawTransport
(
    const word& name,
    const powerLawTransport& ct
)
:
    Thermo(name, ct),
    mu0_(ct.mu0_),
    T0_(ct.T0_),
    k_(ct.k_),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawTransport<Thermo> >
Foam::powerLawTransport<Thermo>::clone() const
{
    return autoPtr<powerLawTransport<Thermo> >
    (
        new powerLawTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawTransport<Thermo> >
Foam::powerLawTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<powerLawTransport<Thermo> >
    (
        new powerLawTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::powerLawTransport<Thermo> >
Foam::powerLawTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<powerLawTransport<Thermo> >
    (
        new powerLawTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::powerLawTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu0_ * pow(T / T0_, k_);
}


template<class Thermo>
inline Foam::scalar Foam::powerLawTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return this->Cp(p, T)*mu(p, T)*rPr_;
}


template<class Thermo>
inline Foam::scalar Foam::powerLawTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return mu(p, T)*rPr_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawTransport<Thermo>& Foam::powerLawTransport<Thermo>::operator=
(
    const powerLawTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu0_ = ct.mu0_;
    T0_ = ct.T0_;
    k_ = ct.k_;
    rPr_ = ct.rPr_;

    return *this;
}


template<class Thermo>
inline void Foam::powerLawTransport<Thermo>::operator+=
(
    const powerLawTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::powerLawTransport<Thermo>::operator-=
(
    const powerLawTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::powerLawTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::powerLawTransport<Thermo> Foam::operator+
(
    const powerLawTransport<Thermo>& ct1,
    const powerLawTransport<Thermo>& ct2
)
{
    notImplemented ("inline Foam::powerLawTransport<Thermo> Foam::operator+");
}


template<class Thermo>
inline Foam::powerLawTransport<Thermo> Foam::operator-
(
    const powerLawTransport<Thermo>& ct1,
    const powerLawTransport<Thermo>& ct2
)
{
    notImplemented ("inline Foam::powerLawTransport<Thermo> Foam::operator-");
}


template<class Thermo>
inline Foam::powerLawTransport<Thermo> Foam::operator*
(
    const scalar s,
    const powerLawTransport<Thermo>& ct
)
{
    notImplemented ("inline Foam::powerLawTransport<Thermo> Foam::operator*");
}


template<class Thermo>
inline Foam::powerLawTransport<Thermo> Foam::operator==
(
    const powerLawTransport<Thermo>& ct1,
    const powerLawTransport<Thermo>& ct2
)
{
    notImplemented ("inline Foam::powerLawTransport<Thermo> Foam::operator==");
}


// ************************************************************************* //
