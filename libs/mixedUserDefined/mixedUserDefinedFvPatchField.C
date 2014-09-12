/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "mixedUserDefinedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
mixedUserDefinedFvPatchField<Type>::mixedUserDefinedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    alfa_(p.size()),
    beta_(p.size()),
    omega0_(p.size())
{}


template<class Type>
mixedUserDefinedFvPatchField<Type>::mixedUserDefinedFvPatchField
(
    const mixedUserDefinedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    alfa_(ptf.alfa_, mapper),
    beta_(ptf.beta_, mapper),
    omega0_(ptf.omega0_, mapper)
{/*
    if (&iF && mapper.hasUnmapped())
    {
        WarningIn
        (
            "mixedUserDefinedFvPatchField<Type>::mixedUserDefinedFvPatchField\n"
            "(\n"
            "    const mixedUserDefinedFvPatchField<Type>&,\n"
            "    const fvPatch&,\n"
            "    const DimensionedField<Type, volMesh>&,\n"
            "    const fvPatchFieldMapper&\n"
            ")\n"
        )   << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
*/
}


template<class Type>
mixedUserDefinedFvPatchField<Type>::mixedUserDefinedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    alfa_("alfa", dict, p.size()),
    beta_("beta", dict, p.size()),
    omega0_("omega0", dict, p.size())
{
    evaluate();
}


template<class Type>
mixedUserDefinedFvPatchField<Type>::mixedUserDefinedFvPatchField
(
    const mixedUserDefinedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    alfa_(ptf.alfa_),
    beta_(ptf.beta_),
    omega0_(ptf.omega0_)
{}


template<class Type>
mixedUserDefinedFvPatchField<Type>::mixedUserDefinedFvPatchField
(
    const mixedUserDefinedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    alfa_(ptf.alfa_),
    beta_(ptf.beta_),
    omega0_(ptf.omega0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void mixedUserDefinedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    alfa_.autoMap(m);
    beta_.autoMap(m);
    omega0_.autoMap(m);
}


template<class Type>
void mixedUserDefinedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const mixedUserDefinedFvPatchField<Type>& mptf =
        refCast<const mixedUserDefinedFvPatchField<Type> >(ptf);

    alfa_.rmap(mptf.alfa_, addr);
    beta_.rmap(mptf.beta_, addr);
    omega0_.rmap(mptf.omega0_, addr);
}


template<class Type>
void mixedUserDefinedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
	Type(pTraits<Type>::one)*( alfa_ / (alfa_+beta_)*omega0_ ) + beta_/(alfa_+beta_)*this->patchInternalField()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > mixedUserDefinedFvPatchField<Type>::snGrad() const
{
    return Type(pTraits<Type>::one)*( alfa_ / (alfa_+beta_)*omega0_ )*this->patch().deltaCoeffs()
		-(alfa_/(alfa_+beta_))*this->patch().deltaCoeffs()*this->patchInternalField() ;
}


template<class Type>
tmp<Field<Type> > mixedUserDefinedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one)*(beta_/(alfa_+beta_));
}


template<class Type>
tmp<Field<Type> > mixedUserDefinedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one)*( alfa_ / (alfa_+beta_)*omega0_ );
}


template<class Type>
tmp<Field<Type> > mixedUserDefinedFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*(alfa_/(alfa_+beta_))*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > mixedUserDefinedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return Type(pTraits<Type>::one)*( alfa_ / (alfa_+beta_)*omega0_ )*this->patch().deltaCoeffs();
}


template<class Type>
void mixedUserDefinedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    alfa_.writeEntry("alfa", os);
    beta_.writeEntry("beta", os);
    omega0_.writeEntry("omega0", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
