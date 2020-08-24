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
    Foam::qgd::constTau

Group
    grpQGDCoeffs

Description
    Class for case constant and uniform \tau distribution

SourceFiles
    constTau.C

\*---------------------------------------------------------------------------*/

#include "QGDCoeffs.H"

#ifndef constTau_H
#define constTau_H

namespace Foam
{

namespace qgd
{

class constTau : public QGDCoeffs
{

    //-
    scalar tau_;

public:

    //-
    TypeName ("constTau");

    //-
    constTau
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~constTau();

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
