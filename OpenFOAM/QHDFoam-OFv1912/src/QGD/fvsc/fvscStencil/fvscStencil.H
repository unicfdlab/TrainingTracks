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
    grpfvscStencil
Class
    Foam::fvsc::fvscStencil
Source file
    fvscStencil.C
\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "regIOobject.H"
#include "runTimeSelectionTables.H"
#include "surfaceFields.H"

#ifndef fvscStencil_H
#define fvscStencil_H

namespace Foam
{

namespace fvsc
{

class fvscStencil : public regIOobject, public refCount
{

protected:
    
    //-
    const fvMesh& mesh_;
    
    //-
    const Time& runTime_;
    
    //-
    surfaceVectorField nf_;
    
    //-
    static PtrList<fvscStencil> stencils_;

public:
 
    //-
    TypeName("fvscStencil");
    
    //-
    declareRunTimeSelectionTable
    (
        autoPtr,
        fvscStencil,
        components,
        (
            const IOobject& io
        ),
        (io)
    );
    
    //- Construct from components
    fvscStencil
    (
        const IOobject& io
    );
    
    //- Return a reference to the selected fvscStencil model
    static autoPtr<fvscStencil> New
    (
        const word& name,
        const fvMesh& mesh
    );
    
    //-
    //static tmp<fvscStencil> lookupOrNew
    static fvscStencil& lookupOrNew
    (
        const word& nname,
        const fvMesh& mesh
    );
    
    //-
    virtual ~fvscStencil();
    
    //-
    virtual tmp<surfaceVectorField> Grad(const volScalarField& vF)
    {
        notImplemented("tmp<surfaceVectorField> Grad(const volScalarField& vF)");
        return tmp<surfaceVectorField>(nullptr);
    }

    //-
    virtual tmp<surfaceTensorField> Grad(const volVectorField& iVF)
    {
        notImplemented("tmp<surfaceTensorField> Grad(const volVectorField& vF)");
        return tmp<surfaceTensorField>(nullptr);
    }

    //-
    virtual tmp<surfaceScalarField> Div(const volVectorField& iVF)
    {
        notImplemented("tmp<surfaceScalarField> Grad(const volVectorField& vF)");
        return tmp<surfaceScalarField>(nullptr);
    }

    //-
    virtual tmp<surfaceVectorField> Div(const volTensorField& iTF)
    {
        notImplemented("tmp<surfaceVectorField> Grad(const volTensorField& vF)");
        return tmp<surfaceVectorField>(nullptr);
    }


    virtual bool writeData(Ostream&) const
    {
        return true;
    }
};

}

}

#endif

//END-OF-FILE

