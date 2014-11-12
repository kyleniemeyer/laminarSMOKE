PtrList<volScalarField> Y;
PtrList<volScalarField> X;
PtrList<volScalarField> C;
Y.resize(thermodynamicsMapXML->NumberOfSpecies());
X.resize(thermodynamicsMapXML->NumberOfSpecies());
C.resize(thermodynamicsMapXML->NumberOfSpecies());

{
	Info<< "Reading species..." << endl;

	for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
	{
		IOobject header
		(
			thermodynamicsMapXML->NamesOfSpecies()[i],
			mesh.time().timeName(),
			mesh,
			IOobject::NO_READ
		);

		// check if field exists and can be read
		if (header.headerOk())
		{
			std::cout << " * " << thermodynamicsMapXML->NamesOfSpecies()[i] << " -> from file..." << std::endl;
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
						IOobject::NO_WRITE
					),
					mesh
				)
			);
		}
		else
		{
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
						IOobject::NO_WRITE
					),
					Ydefault
				)
			);
		}

		X.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"x::" + thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				Y[i]
			)
		);

		C.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					"c::" + thermodynamicsMapXML->NamesOfSpecies()[i],
					mesh.time().timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionSet(0, -3, 0, 0, 1)
			)
		);
	}
}
