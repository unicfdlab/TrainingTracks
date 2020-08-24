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
Group grpfvsc
    This group contains common part of QGD solvers.
Class
    Foam::fvsc
Description 
    Methods calculating of differential operators
\*---------------------------------------------------------------------------*/


#include "fvsc.H"
#include "fvscStencil.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "prismMatcher.H"

namespace Foam
{
namespace fvsc
{
    word fvscOpName(const Foam::fvMesh& mesh, Foam::word termName);
}
}

Foam::word
Foam::fvsc::fvscOpName(const Foam::fvMesh& mesh, Foam::word termName)
{
    word opname = "none";
    if (mesh.schemesDict().subDict("fvsc").found(termName))
    {
        mesh.schemesDict().subDict("fvsc").lookup(termName) >> opname;
    }
    else
    {
        mesh.schemesDict().subDict("fvsc").lookup("default") >> opname;
    }
     
    if (((opname == "leastSquares") or (opname == "leastSquaresOpt")) and (mesh.nGeometricD() == 3))
    {
	FatalErrorIn("Foam::fvsc::fvscOpName") << "Can't use leastSquares or leastSquaresOpt in 3D case." << nl << exit(FatalError);
    }
    
    if (opname == "GaussVolPoint")
    {
        const labelList wedgeBC = mesh.boundaryMesh().findIndices("wedge", true);
        if (wedgeBC.size() > 0)
        {
            prismMatcher prism;
            forAll (mesh.cells(), cellI)
            {
                if (prism.isA(mesh, cellI))
                {
                    FatalErrorInFunction
                        << "GaussVolPoint scheme does not support solving axisymmetric cases with wedge BC and prism cells." << endl
                        << "Try to set leastSquares scheme."
                        << exit(FatalError);
                }
            }
        }
    }

    return opname;
}

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::grad(const volScalarField& vf)
{
    word tname = "grad(" + vf.name() + ")";
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        fvscOpName(vf.mesh(),tname),
        vf.mesh()
    );
    
    tmp<surfaceVectorField> tGrad(Stencil.Grad(vf));
    
    return tGrad;
}

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::grad(const tmp<volScalarField>& tvf)
{
    return Foam::fvsc::grad(tvf());
}

Foam::tmp<Foam::surfaceTensorField>
Foam::fvsc::grad(const volVectorField& vf)
{
    word tname = "grad(" + vf.name() + ")";
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        fvscOpName(vf.mesh(),tname),
        vf.mesh()
    );
    
    return Stencil.Grad(vf);
}

Foam::tmp<Foam::surfaceTensorField>
Foam::fvsc::grad(const tmp<volVectorField>& tvf)
{
    return Foam::fvsc::grad(tvf());
}

Foam::tmp<Foam::surfaceScalarField>
Foam::fvsc::div(const volVectorField& vf)
{
    word tname = "div(" + vf.name() + ")";
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        fvscOpName(vf.mesh(),tname),
        vf.mesh()
    );
    
    return Stencil.Div(vf);
}

Foam::tmp<Foam::surfaceScalarField>
Foam::fvsc::div(const tmp<volVectorField>& tvf)
{
    return Foam::fvsc::div(tvf());
}

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::div(const volTensorField& vf)
{
    word tname = "div(" + vf.name() + ")";
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        fvscOpName(vf.mesh(),tname),
        vf.mesh()
    );
    
    return Stencil.Div(vf);
}

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::div(const tmp<volTensorField>& tvf)
{
    return Foam::fvsc::div(tvf());
}


//END-OF-FILE


