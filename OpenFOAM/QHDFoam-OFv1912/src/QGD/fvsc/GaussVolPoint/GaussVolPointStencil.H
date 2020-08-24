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
Group 
    grpGaussVolPoint
Class
    Foam::fvsc::GaussVolPoint::GaussVolPointStencil
SourceFile
    GaussVolPointStencil.C 
\*---------------------------------------------------------------------------*/

#ifndef GaussVolPointStencil_H
#define GaussVolPointStencil_H

#include "fvscStencil.H"
#include "volFields.H"
#include "GaussVolPointBase.H"

namespace Foam
{

namespace fvsc
{

/*---------------------------------------------------------------------------*\
                    Class reduced Declaration
\*---------------------------------------------------------------------------*/

class GaussVolPoint
:
    public fvscStencil,
    public GaussVolPointBase
{

public:

        TypeName ("GaussVolPoint");

    // Constructors

        //- Construct from IOobject. 
        GaussVolPoint(const IOobject&);

        //- Destructor
        ~GaussVolPoint();
        
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

