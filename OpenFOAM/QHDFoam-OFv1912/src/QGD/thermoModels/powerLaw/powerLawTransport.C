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
    
Group
    grpPowerLawQGD

\*---------------------------------------------------------------------------*/

#include "powerLawTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::powerLawTransport<Thermo>::powerLawTransport(Istream& is)
:
    Thermo(is),
    mu0_(readScalar(is)),
    T0_(readScalar(is)),
    k_(readScalar(is)),
    rPr_(1.0/readScalar(is))
{
    is.check("powerLawTransport::powerLawTransport(Istream& is)");
}


template<class Thermo>
Foam::powerLawTransport<Thermo>::powerLawTransport(const dictionary& dict)
:
    Thermo(dict),
    mu0_(readScalar(dict.subDict("transport").lookup("mu0"))),
    T0_(readScalar(dict.subDict("transport").lookup("T0"))),
    k_(readScalar(dict.subDict("transport").lookup("k"))),
    rPr_(1.0/readScalar(dict.subDict("transport").lookup("Pr")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::powerLawTransport<Thermo>::powerLawTransport::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("mu0", mu0_);
    dict.add("T0", T0_);
    dict.add("k", k_);
    dict.add("Pr", 1.0/rPr_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const powerLawTransport<Thermo>& ct)
{
    operator<<(os, static_cast<const Thermo&>(ct));
    os<< tab << ct.mu0_ << tab << ct.T0_
    << tab << ct.k_ << tab << 1.0/ct.rPr_;

    os.check("Ostream& operator<<(Ostream&, const powerLawTransport&)");

    return os;
}


// ************************************************************************* //
