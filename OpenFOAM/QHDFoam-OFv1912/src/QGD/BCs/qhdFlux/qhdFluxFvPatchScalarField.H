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
    Foam::qhdFluxFvPatchScalarField

Group
    grpQGDBCs grpWallBoundaryConditions

Description
    This boundary condition sets the pressure gradient to the provided value
    such that the flux on the boundary is that specified by the velocity
    boundary condition.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            qhdFlux;
    }
    \endverbatim

SourceFiles
    qhdFluxFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef qhdFluxFvPatchScalarField_H
#define qhdFluxFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class qhdFluxFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class qhdFluxFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{


public:

    //- Runtime type information
    TypeName("qhdFlux");


    // Constructors

        //- Construct from patch and internal field
        qhdFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        qhdFluxFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given qhdFluxFvPatchScalarField onto
        //  a new patch
        qhdFluxFvPatchScalarField
        (
            const qhdFluxFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        qhdFluxFvPatchScalarField
        (
            const qhdFluxFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new qhdFluxFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        qhdFluxFvPatchScalarField
        (
            const qhdFluxFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new qhdFluxFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        ////- Update the patch pressure gradient field from the given snGradp
        //virtual void updateSnGrad(const scalarField& snGradp);

        //- Update the patch pressure gradient field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "volFields.H"
//
//namespace Foam
//{
//    template<class GradBC>
//    inline void setSnGrad
//    (
//        volScalarField::Boundary& bf,
//        const FieldField<fvsPatchField, scalar>& snGrad
//    )
//    {
//        forAll(bf, patchi)
//        {
//            if (isA<GradBC>(bf[patchi]))
//            {
//                refCast<GradBC>(bf[patchi]).updateSnGrad(snGrad[patchi]);
//            }
//        }
//    }
//
//    template<class GradBC>
//    inline void setSnGrad
//    (
//        volScalarField::Boundary& bf,
//        const tmp<FieldField<fvsPatchField, scalar>>& tsnGrad
//    )
//    {
//        setSnGrad<GradBC>(bf, tsnGrad());
//    }
//}
//

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
