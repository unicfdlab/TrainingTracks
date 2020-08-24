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
Description
    This is a method for calculation the differential operators without 
    tangential derivatives. They are further used in the calculation of 
    the QGD terms.
\*---------------------------------------------------------------------------*/

#include "fvscStencil.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledFvsPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvsc
{
    defineTypeNameAndDebug(fvscStencil, 0);
    defineRunTimeSelectionTable(fvscStencil, components);
}
}

namespace Foam
{
namespace fvsc
{

PtrList<fvscStencil> fvscStencil::stencils_(0);

autoPtr<fvscStencil> fvscStencil::New
(
    const word& fvscType,
    const fvMesh& mesh
)
{
    Info<< "Selecting finite volume surface calculus stencil type " << fvscType << endl;
    
    componentsConstructorTable::iterator cstrIter =
        componentsConstructorTablePtr_->find(fvscType);
    
    if (cstrIter == componentsConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "fvscStencil::New(const word&, const fvMesh&)"
        )   << "Unknown Model type " << fvscType << nl << nl
        << "Valid model types are:" << nl
        << componentsConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }
    
    return autoPtr<fvscStencil>
    (
        cstrIter()
        (
            IOobject
            (
                fvscType,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        )
    );
}

//tmp<fvscStencil> fvscStencil::lookupOrNew
fvscStencil& fvscStencil::lookupOrNew
(
    const word& name,
    const fvMesh& mesh
)
{
    if (!mesh.thisDb().foundObject<fvscStencil>(name))
    {
        stencils_.append
        (
            fvscStencil::New(name,mesh)
        );
        stencils_.last().checkIn();
    }
    
    return
        const_cast<fvscStencil&>
        (
            mesh.thisDb().lookupObject<fvscStencil>(name)
        );
}
//
fvscStencil::fvscStencil(const IOobject& io)
:
    regIOobject(io, false),
    refCount(),
    mesh_(refCast<const fvMesh>(io.db())),
    runTime_(mesh_.time()),
    nf_
    (
        mesh_.Sf() / mesh_.magSf()
    )
{
    nf_.rename("nf");
}

fvscStencil::~fvscStencil()
{
}

}; //namespace fvsc

}; //namespace Foam


//END-OF-FILE

