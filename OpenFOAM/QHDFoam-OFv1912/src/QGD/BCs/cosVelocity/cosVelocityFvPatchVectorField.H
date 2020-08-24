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
    Foam::cosVelocityFvPatchVectorField

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
        type            cosVelocity;
    }
    \endverbatim

SourceFiles
    cosVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cosVelocityFvPatchVectorField_H
#define cosVelocityFvPatchVectorField_H

#include "mixedFvPatchFields.H"
#include "coordinateSystem.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class advectiveFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cosVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{
protected:

    // Private data

        //- 
        vector A_;
        
        //-
        scalar H_;
        
        //-
        vector Hdirection_;

        //- 
        scalar omega0_;
        
        //-
        scalar phi0_;
        
        //-
        scalar minZ_;

public:

    //- Runtime type information
    TypeName("cosVelocity");


    // Constructors

        //- Construct from patch and internal field
        cosVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cosVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given cosVelocityFvPatchVectorField
        //  onto a new patch
        cosVelocityFvPatchVectorField
        (
            const cosVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cosVelocityFvPatchVectorField
        (
            const cosVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField > clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new cosVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cosVelocityFvPatchVectorField
        (
            const cosVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new cosVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access
        void setPatchVelocities(vectorField& pV);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
        
        //-
        virtual void autoMap (const fvPatchFieldMapper& );
        
        //-
        virtual void rmap (const fvPatchVectorField& pf, const labelList& ll);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //
