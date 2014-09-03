{
    if(energyEquation == true)
    {
		sumCpSpecies *= 0.;
		if (iMassDiffusionInEnergyEquation == true)
		{

			for (label i=0; i<Y.size(); i++)
			{
				volScalarField& Yi 		= Y[i];			
				volScalarField& Dmixi 	= Dmix[i];
				volScalarField& Cpi 	= CpSpecies[i];
		
				sumCpSpecies += rho*Cpi*Dmixi*(fvc::grad(Yi)&fvc::grad(T));
			}
		}

		fvScalarMatrix TEqn
		(
			  cp*fvm::ddt(rho,T)
			+ cp*fvm::div(phi, T)		    
		 ==
		    fvm::laplacian(lambda,T) + 
		    sumCpSpecies +
			cp*fvOptions(rho, T)
		);

		TEqn.relax();
		fvOptions.constrain(TEqn);
		TEqn.solve(mesh.solver("T"));
		fvOptions.correct(T);	
    }

    Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;

	// Spark
	spark.SetSpark(runTime.value(), T);
}
