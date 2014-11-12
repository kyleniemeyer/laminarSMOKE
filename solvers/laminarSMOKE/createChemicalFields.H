
Info << endl;
Info<< " * Reading chemical species mole fractions or mass fractions..." << endl;

PtrList<volScalarField> Y(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> Dmix(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> HSpecies(thermodynamicsMapXML->NumberOfSpecies());
PtrList<volScalarField> FormationRates(outputFormationRatesIndices.size());

#if STEADYSTATE == 1
PtrList<volScalarField> sourceImplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
PtrList<volScalarField> sourceExplicit(thermodynamicsMapXML->NumberOfSpecies()+1);
#endif

// Loop over all the species in the kinetic mechanism
for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
{
	Info << "   Reading species: " << thermodynamicsMapXML->NamesOfSpecies()[i] << " > ";

	IOobject header
	(
		thermodynamicsMapXML->NamesOfSpecies()[i],
		mesh.time().timeName(),
		mesh,
		IOobject::NO_READ
	);

	// Sets of diffusivities and specific heats of species
	{
		Dmix.set
		(
			i,
			new volScalarField
			(
			IOobject
			(
				"gas::Dmix_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh,
			dimensionSet(0, 2, -1, 0, 0)
			)
		);

		HSpecies.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
				"thermo::H_"+ thermodynamicsMapXML->NamesOfSpecies()[i],
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
				),
				mesh,
				dimensionSet(0, 2, -2, 0, 0)
			)
		);
	}

	// check if field exists and can be read
	if (header.headerOk())
	{
		Info << thermodynamicsMapXML->NamesOfSpecies()[i] << endl;
		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);
	}
	else
	{
		Info << "Ydefault" << endl;

		volScalarField Ydefault
		(
			IOobject
			(
				"Ydefault",
				mesh.time().timeName(),
				mesh,
				IOobject::MUST_READ,
				IOobject::NO_WRITE
			),
			mesh
		);

		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				Ydefault
			)
		);
	}

	if (outputFormationRatesIndices.size() != 0)
	{
		for (int i=0;i<outputFormationRatesIndices.size();i++)
		{
			FormationRates.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"R::" + thermodynamicsMapXML->NamesOfSpecies()[outputFormationRatesIndices(i)],
						mesh.time().timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					mesh,
					dimensionedScalar("R", dimensionSet(1, -3, -1, 0, 0), 0.0)
				)
			);
		}
	}

	#if STEADYSTATE == 1

	for (int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
	{
		sourceImplicit.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"SI::" + thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("SI", dimensionSet(1, -3, -1, 0, 0), 0.0)
			)
		);

		sourceExplicit.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"SE::" + thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionedScalar("SE", dimensionSet(1, -3, -1, 0, 0), 0.0)
			)
		);
	}

	sourceImplicit.set
	(
		thermodynamicsMapXML->NumberOfSpecies(),
		new volScalarField
		(
			IOobject
			(
				"SI::T",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("SIT", dimensionSet(1, -1, -3, -1, 0), 0.0)
		)
	);

	sourceExplicit.set
	(
		thermodynamicsMapXML->NumberOfSpecies(),
		new volScalarField
		(
			IOobject
			(
				"SE::T",
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			mesh,
			dimensionedScalar("SET", dimensionSet(1, -1, -3, 0, 0), 0.0)
		)
	);

	#endif
}
 

PtrList<volScalarField> Dsoret(transportMapXML->iThermalDiffusionRatios().size());
if (soretEffect == true)
{
	for(unsigned int i=0;i<transportMapXML->iThermalDiffusionRatios().size();i++)
	{ 
		unsigned int index = transportMapXML->iThermalDiffusionRatios()[i];
		indexSoret[index-1] = i;
		soretEffectList[index-1] = true;

		std::cout << "Soret effect: species " 	<< i+1 << "/" << transportMapXML->iThermalDiffusionRatios().size() 
							<< " " << thermodynamicsMapXML->NamesOfSpecies()[index-1] << " (" << index << ")" << std::endl;

		Dsoret.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"gas::Dsoret_"+ thermodynamicsMapXML->NamesOfSpecies()[index-1],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::AUTO_WRITE
				),
				mesh,
				dimensionSet(0, 2, -1, 0, 0)
			)
		);
	}
}

#include "createFluxes.H"
