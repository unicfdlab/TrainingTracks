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
    Foam::qgd::H2bynuQHD

Group
    grpQGDCoeffs

Description
    Class for one of possible ways of tau calculating.

$$    \tau_QGD = \alpha \cdot \frac{h^2}{\nu}  $$

    where $\tau_QGD$ - regulatisation parameter, $\alpha$ - parameter
    which depends by user in constant foulder, $h$ - size of mesh,
    $\nu$ - kinetic viscosity.

SourceFiles
    H2bynuQHD.C

\*---------------------------------------------------------------------------*/

#include "QGDCoeffs.H"

#ifndef H2bynuQHD_H
#define H2bynuQHD_H

namespace Foam
{

namespace qgd
{

class H2bynuQHD : public QGDCoeffs
{

public:

    //-
    TypeName ("H2bynuQHD");

    //-
    H2bynuQHD
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~H2bynuQHD();

    //-
    void correct(const QGDThermo&);
};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//
