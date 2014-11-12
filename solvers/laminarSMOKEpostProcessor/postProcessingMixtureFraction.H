{
	//- Mixture fraction
	Info<< "Creating mixture fraction\n" << endl;
	volScalarField Z_
	(
	   	IOobject
	   	(
	     	   	"post_Z",
	     	   	runTime.timeName(),
	    	   	mesh,
	    	    	IOobject::NO_READ,
	    	    	IOobject::NO_WRITE
	   	),
	   	 mesh,
	  	  dimensionSet(0, 0, 0, 0, 0)
	);

	// Calculations
	{
		std::ofstream fTable("Table.out", std::ios::out);
		fTable.setf(std::ios::scientific);

		std::vector<std::string> fuel_names(1); fuel_names[0] = "CH4";
		std::vector<std::string> oxidizer_names(2); oxidizer_names[0] = "O2"; oxidizer_names[1] = "N2";
		std::vector<double> mass_fuel(1); mass_fuel[0] = 1.;
		std::vector<double> mass_oxidizer(2); mass_oxidizer[0] = 0.232; mass_oxidizer[1] = 0.768;

		label ns = Y_.size();

		OpenSMOKE::OpenSMOKEVectorDouble y(ns);
		OpenSMOKE::OpenSMOKEVectorDouble x(ns);
		std::vector<double> y_stl(ns);

		const scalarField& pCells = p.internalField();
		const scalarField& TCells = T.internalField();
		const scalarField& ZCells = Z_.internalField();

		forAll(TCells, celli)
		{	
			double mw;

			// Extract the mean mass fractions
			for(unsigned int i=0;i<ns;i++)
				y[i+1] = Y_[i].internalField()[celli];

			// Convert to mole fractions
			thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

			MW_.internalField()[celli] = mw;
			for(unsigned int i=0;i<ns;i++)
				X_[i].internalField()[celli] = x[i+1];

			for(unsigned int i=0;i<ns;i++)
				y_stl[i] = y[i+1];
			Z_.internalField()[celli] = thermodynamicsMapXML->GetMixtureFractionFromMassFractions(y_stl, fuel_names, mass_fuel,oxidizer_names, mass_oxidizer);

			if (Z_.internalField()[celli]>=0. && Z_.internalField()[celli]<=1.)
			{
				fTable << T.internalField()[celli] << " " << Z_.internalField()[celli] << " ";
				for(unsigned int i=0;i<ns;i++)
					fTable << y[i+1] << " ";
				fTable << std::endl;
			}
		}

		forAll(T.boundaryField(), patchi)
		{
			fvPatchScalarField& pp = p.boundaryField()[patchi];
			fvPatchScalarField& pT = T.boundaryField()[patchi];
			fvPatchScalarField& pMW = MW_.boundaryField()[patchi];
			fvPatchScalarField& pZ = Z_.boundaryField()[patchi];
	
			forAll(pT, facei)
			{
				double mw;

				// Sets the maps
				thermodynamicsMapXML->SetTemperature(pT[facei]);
				thermodynamicsMapXML->SetPressure(pp[facei]);

				// Extract the mass fractions
				for(unsigned int i=1;i<=ns;i++)
					y[i] = Y_[i-1].boundaryField()[patchi][facei];

				// Convert to mole fractions
				thermodynamicsMapXML->MoleFractions_From_MassFractions(x, mw, y);

				pMW[facei] = mw;
				for(unsigned int i=1;i<=ns;i++)
					X_[i-1].boundaryField()[patchi][facei] = x[i];

				for(unsigned int i=0;i<ns;i++)
					y_stl[i] = y[i+1];
				pZ[facei] = thermodynamicsMapXML->GetMixtureFractionFromMassFractions(y_stl, fuel_names, mass_fuel,oxidizer_names, mass_oxidizer);

				if (pZ[facei]>=0. && pZ[facei]<=1.)
				{
					fTable << pT[facei] << " " << pZ[facei] << " ";
					for(unsigned int i=0;i<ns;i++)
						fTable << y[i+1] << " ";
					fTable << std::endl;
				}
			}
		}

		fTable.close();
		MW_.write();
		Z_.write();			
	}
}
