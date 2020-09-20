/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::RNSfixedRhoFvPatchScalarField

Description
    Foam::RNSfixedRhoFvPatchScalarField

    This boundary condition provides a fixed density inlet condition for
    compressible solvers, where the density of calculated using:

    \f[
        \rho = \psi p
    \f]

    where
    \vartable
        p     | pressure [Pa]
        \rho  | density [kg/m3]
    \endvartable


Usage
    \table
        Property     | Description                | Required    | Default value
        p            | Pressure field name        | no          | p
        psi          | Compressibility field name | no          | thermo:psi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedRho;
    }
    \endverbatim

SourceFiles
    RNSfixedRhoFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef RNSfixedRhoFvPatchScalarField_H
#define RNSfixedRhoFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class RNSfixedRhoFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class RNSfixedRhoFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

private:

    // Private data

        //- Pressure field name, default = "p"
        word pName_;

        //- Compressibility field name, default = "thermo:psi"
        word psiName_;


public:

    //- Runtime type information
    TypeName("RNSfixedRho");


    // Constructors

        //- Construct from patch and internal field
        RNSfixedRhoFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        RNSfixedRhoFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given RNSfixedRhoFvPatchScalarField
        // onto a new patch
        RNSfixedRhoFvPatchScalarField
        (
            const RNSfixedRhoFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        RNSfixedRhoFvPatchScalarField
        (
            const RNSfixedRhoFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new RNSfixedRhoFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        RNSfixedRhoFvPatchScalarField
        (
            const RNSfixedRhoFvPatchScalarField&,
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
                new RNSfixedRhoFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
