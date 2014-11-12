{
    if(energyEquation == true)
    {
		massDiffusionInEnergyEquation *= 0.;
		if (iMassDiffusionInEnergyEquation == true)
		{
			for (label i=0; i<Y.size(); i++)
				massDiffusionInEnergyEquation -= fvc::div(J[i],HSpecies[i], "div(Ji,hsi)") - 
								 HSpecies[i] * fvc::div(J[i]);
		}

		if (radiativeHeatTransfer.radiationModel() == RADIATION_NONE)
		{
			fvScalarMatrix TEqn
			(
				cp*fvm::div(phi, T)		    
				 ==
		    		fvm::laplacian(lambda,T) + 
		   	 	massDiffusionInEnergyEquation +
                                sourceExplicit[Y.size()] + fvm::Sp(sourceImplicit[Y.size()],T) +
		    		cp*fvOptions(rho, T)
			);

			TEqn.relax();
			fvOptions.constrain(TEqn);
			TEqn.solve(mesh.solver("T"));
			fvOptions.correct(T);	
		}
		else
		{
			fvScalarMatrix TEqn
			(
				cp*fvm::div(phi, T)		    
				 ==
		    		fvm::laplacian(lambda,T) + 
		   	 	massDiffusionInEnergyEquation - 
                                Qrad +
                                sourceExplicit[Y.size()] + fvm::Sp(sourceImplicit[Y.size()],T) +
		    		cp*fvOptions(rho, T)
			);

			TEqn.relax();
			fvOptions.constrain(TEqn);
			TEqn.solve(mesh.solver("T"));
			fvOptions.correct(T);
		}
    }

    Info<< " * T gas min/max (after transport)   = " << min(T).value() << ", " << max(T).value() << endl;

    // Spark
    spark.SetSpark(runTime.value(), T);
}
