tmp<volScalarField> Rho(volScalarField& p, volScalarField& psi)
{
	return p*psi;
}

void normalizeMassFractions(Eigen::VectorXd& omega_plus_temperature, const label celli, const double massFractionsTol)
{
	double sumFractions = 0.;
    for(int i=0; i < omega_plus_temperature.size()-1; i++)
            sumFractions += omega_plus_temperature(i);

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of mass-fractions = " << sumFractions << " in cell " << celli <<endl;

        //for(int i=0; i < omega_plus_temperature.size()-1; i++)
        //    cout << i << ")\t" << omega_plus_temperature(i) << endl;
        
        if (time == 0)
        {
            Info << "\nFATAL ERROR: sum of Yi is not 1" <<endl;
            Info << "Check internal field on cell: " << celli <<endl;
            abort();
        }    
    }
	
	 for(int i=0; i < omega_plus_temperature.size()-1; i++)
            omega_plus_temperature(i) /= sumFractions;
}

void normalizeMassFractions(OpenSMOKE::OpenSMOKEVectorDouble& omega_plus_temperature, const label celli, const double massFractionsTol)
{
	double sumFractions = 0.;
    for(int i=1; i <= omega_plus_temperature.Size()-1; i++)
            sumFractions += omega_plus_temperature[i];

    if (sumFractions > 1.+massFractionsTol || sumFractions < 1.-massFractionsTol)
    {
        Info << "WARNING: sum of mass-fractions = " << sumFractions << " in cell " << celli <<endl;

        //for(int i=1; i <= omega_plus_temperature.Size()-1; i++)
        //    cout << i << ")\t" << omega_plus_temperature[i] << endl;
        
        if (time == 0)
        {
            Info << "\nFATAL ERROR: sum of Yi is not 1" <<endl;
            Info << "Check internal field on cell: " << celli <<endl;
            abort();
        }    
    }
	
	 for(int i=1; i <= omega_plus_temperature.Size()-1; i++)
            omega_plus_temperature[i] /= sumFractions;
}
