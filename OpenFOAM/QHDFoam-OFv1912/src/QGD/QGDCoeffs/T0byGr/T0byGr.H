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
    Foam::qgd::T0byGr

Group
    grpQGDCoeffs

Description
    Class for to calculate $\tau$ as inverse proportional of 
    Grashoff number.
    in this case $\tau$ parameter is calculated by

    $$ \tau = \frac{T_0}{Gr} $$

    Where:

    $\tau$ is regularisation parameter, $T_0$ is some
    characteristic time.
    
    $Gr$ is Grashof number:

    $$\frac{g \beta \Delta TH}{\nu^2}$$

    Where:

    $\beta$ is the coefficient of thermal expansion
    $\Delta T$ is deiffetence of temperatures / salinities
    $\nu$ is the kinematic viscosity

SourceFiles
    T0byGr.C

\*---------------------------------------------------------------------------*/

#include "QGDCoeffs.H"

#ifndef T0byGr_H
#define T0byGr_H

namespace Foam
{

namespace qgd
{

class T0byGr : public QGDCoeffs
{

    //-
    scalar Gr_;
    
    //-
    scalar T0_;

public:

    //-
    TypeName ("T0byGr");

    //-
    T0byGr
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~T0byGr();

    //-
    void correct(const QGDThermo&);

private:

    //-

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//
