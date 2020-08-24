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
    Foam::fvsc::GaussVolPoint::GaussVolPoint1D
Description
    This is a method for approximating derivatives of tangents to a face (1D case). 
    They are further used in the calculation of the QGD terms.
\*---------------------------------------------------------------------------*/
#include "GaussVolPointBase1D.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"

Foam::fvsc::GaussVolPointBase1D::GaussVolPointBase1D(const fvMesh& mesh)
:
    nfRef_(mesh.thisDb().lookupObject<surfaceVectorField>("nf"))
{
};


Foam::fvsc::GaussVolPointBase1D::~GaussVolPointBase1D()
{
}

void Foam::fvsc::GaussVolPointBase1D::faceGrad(const volScalarField& f, surfaceVectorField& gradf)
{
    if (f.mesh().nGeometricD() == 1)
    {
        gradf = nfRef_() * fvc::snGrad(f);
    }
};

void Foam::fvsc::GaussVolPointBase1D::faceGrad(const volVectorField& f, surfaceTensorField& gradf)
{
    if (f.mesh().nGeometricD() == 1)
    {
        gradf = nfRef_() * fvc::snGrad(f);
    }
};

void Foam::fvsc::GaussVolPointBase1D::faceDiv(const volVectorField& f, surfaceScalarField& divf)
{
    if (f.mesh().nGeometricD() == 1)
    {
        divf = nfRef_() & fvc::snGrad(f);
    }
};

void Foam::fvsc::GaussVolPointBase1D::faceDiv(const volTensorField& f, surfaceVectorField& divf)
{
    if (f.mesh().nGeometricD() == 1)
    {
        divf = nfRef_() & fvc::snGrad(f);
    }
}

//
//END-OF-FILE
//


