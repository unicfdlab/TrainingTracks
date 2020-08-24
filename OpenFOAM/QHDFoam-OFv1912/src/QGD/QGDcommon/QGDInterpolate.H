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

Description
    Creates interpolation instances templated for QGD solver

SourceFiles
    QGDInterpolate.H

\*---------------------------------------------------------------------------*/
namespace Foam
{

template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdInterpolate(const GeometricField<T, Foam::fvPatchField, Foam::volMesh > &psi)
{
    if
    (
        psi.mesh().schemesDict().subDict("interpolationSchemes").found
        ("interpolate("+psi.name()+")")
    )
    {
        return fvc::interpolate(psi);
    }
    else if 
    (
        psi.mesh().schemesDict().subDict("interpolationSchemes").found
        ("default")
    )
    {
        if
        (
            psi.mesh().schemesDict().subDict("interpolationSchemes").template get<word>
            ("default") == "none"
        )
        {
            return linearInterpolate(psi);
        }
        return fvc::interpolate(psi);
    }
    return linearInterpolate(psi);
}

template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdInterpolate(const tmp<GeometricField<T, Foam::fvPatchField, Foam::volMesh >> &tpsi)
{
    return qgdInterpolate(tpsi());
}

template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdFlux
(
    const GeometricField<scalar, Foam::fvsPatchField, Foam::surfaceMesh>& flux,
    const GeometricField<T, Foam::fvPatchField, Foam::volMesh > &psi,
    const GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh>& psif,
    const word fluxName
)
{
    if (psi.mesh().schemesDict().subDict("divSchemes").found(fluxName))
    {
        //if
        //(
        //    psi.mesh().schemesDict().subDict("divSchemes").template get<word>(fluxName)
        //    ==
        //    "linear"
        //)
        //{
        //    return flux*psif;
        //}
        return fvc::flux
        (
            flux,
            psi,
            fluxName
        );
    }
    return flux*psif;
}

template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdFlux
(
    const GeometricField<scalar, Foam::fvsPatchField, Foam::surfaceMesh>& flux,
    const GeometricField<T, Foam::fvPatchField, Foam::volMesh > &psi,
    const GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh>& psif
)
{
    word fluxName = "div("+flux.name()+","+psi.name()+")";
    return qgdFlux(flux,psi,psif,fluxName);
}


template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdFlux
(
    const tmp<GeometricField<scalar, Foam::fvsPatchField, Foam::surfaceMesh>>& flux,
    const tmp<GeometricField<T, Foam::fvPatchField, Foam::volMesh >> &psi,
    const tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh>>& psif
)
{
    return qgdFlux(flux(),psi(),psif());
}

}

//
//END-OF-FILE
//

