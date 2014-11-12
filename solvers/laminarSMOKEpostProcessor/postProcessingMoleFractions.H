// Post processing: mole fractions
{
	label ns = Y.size();

	OpenSMOKE::OpenSMOKEVectorDouble y(ns);
	OpenSMOKE::OpenSMOKEVectorDouble x(ns);

	const scalarField& TCells = T.internalField();
	
	forAll(TCells, celli)
	{	
		double mw;

		// Extract the mean mass fractions
		for(unsigned int i=0;i<ns;i++)
			y[i+1] = Y[i].internalField()[celli];

		// Convert to mole fractions
		thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

		for(unsigned int i=0;i<ns;i++)
			X[i].internalField()[celli] = x[i+1];
	}

	forAll(T.boundaryField(), patchi)
	{
		fvPatchScalarField& pT = T.boundaryField()[patchi];

		forAll(pT, facei)
		{
			double mw;

			// Extract the mass fractions
			for(unsigned int i=1;i<=ns;i++)
				y[i] = Y[i-1].boundaryField()[patchi][facei];

			// Convert to mole fractions
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

			for(unsigned int i=1;i<=ns;i++)
				X[i-1].boundaryField()[patchi][facei] = x[i];
		}
	}

	for (unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		X[i].write();		
}
