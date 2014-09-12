/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "backwardDiffusionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedUserDefinedFvPatchScalarField(p, iF)
{
	alfa() = 0.0;
	beta() = 0.0;
	omega0() = 0.0;
}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const backwardDiffusionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedUserDefinedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedUserDefinedFvPatchScalarField(p, iF)
{
    // Set the nominal value
    omega0() = scalarField("omega0", dict, p.size());

    // Fixed value condition is forced
    alfa() = 1000.;
    const double Dmix = 1e-10;
    beta() = Dmix*this->patch().deltaCoeffs();

    // Read value if available
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate();
    }
}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const backwardDiffusionFvPatchScalarField& tppsf
)
:
    mixedUserDefinedFvPatchScalarField(tppsf)
{}


Foam::backwardDiffusionFvPatchScalarField::backwardDiffusionFvPatchScalarField
(
    const backwardDiffusionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedUserDefinedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::backwardDiffusionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void Foam::backwardDiffusionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedUserDefinedFvPatchScalarField::rmap(ptf, addr);

//    const backwardDiffusionFvPatchScalarField& tiptf =
//        refCast<const backwardDiffusionFvPatchScalarField>(ptf);
}


void Foam::backwardDiffusionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const volVectorField& U = db().lookupObject<volVectorField>("U");
    tmp<vectorField> n = patch().nf();
    alfa() = -(n & U.boundaryField()[patchi]);

    nameInternal_ = dimensionedInternalField().name();
    const volScalarField& Dmix = db().lookupObject<volScalarField>("gas::Dmix_" + nameInternal_);
    beta() = Dmix.boundaryField()[patchi]*this->patch().deltaCoeffs();

    if (debug)
    {
    }

    mixedUserDefinedFvPatchScalarField::updateCoeffs();
}

void Foam::backwardDiffusionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    omega0().writeEntry("omega0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        backwardDiffusionFvPatchScalarField
    );
}

// ************************************************************************* //
