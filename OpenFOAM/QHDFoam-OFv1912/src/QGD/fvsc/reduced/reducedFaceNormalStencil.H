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
Group grpreduced
    This group contains common part of QGD solvers.
Class
    Foam::fvsc::reduced::reducedFaceNormalStencil
Description 
    Methods calculating of differential operators without tangential direvatives
Source files
    reducedFaceNormalStencil.C   
\*---------------------------------------------------------------------------*/

#ifndef reducedFaceNormalStencil_H
#define reducedFaceNormalStencil_H

#include "fvscStencil.H"
#include "volFields.H"

namespace Foam
{

namespace fvsc
{

/*---------------------------------------------------------------------------*\
                    Class reduced Declaration
\*---------------------------------------------------------------------------*/

class reduced
:
    public fvscStencil
{

public:

        TypeName ("reduced");

    // Constructors

        //- Construct from IOobject. 
        reduced(const IOobject&);

        //- Destructor
        ~reduced();
        
        //-
        tmp<surfaceVectorField> Grad(const volScalarField& iF);
        
        //-
        tmp<surfaceTensorField> Grad(const volVectorField& iVF);
        
        //-
        tmp<surfaceScalarField> Div(const volVectorField& iVF);
        
        //-
        tmp<surfaceVectorField> Div(const volTensorField& iTF);
};

} //namespace fvsc

} //namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

//END-OF-FILE

