// 
void reactionSourceTerms(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap_,
				const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0, OpenSMOKE::OpenSMOKEVectorDouble& S)
{
	const unsigned int NC_ = thermodynamicsMap_.NumberOfSpecies();
	OpenSMOKE::OpenSMOKEVectorDouble omega_(NC_);
	OpenSMOKE::OpenSMOKEVectorDouble x_(NC_);
	OpenSMOKE::OpenSMOKEVectorDouble R_(NC_);
	OpenSMOKE::OpenSMOKEVectorDouble c_(NC_);
	
	for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = max(y[i], 0.);
	const double T = y[NC_+1];

	// Calculates the pressure and the concentrations of species
	double MW_ = 0.;
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
	const double cTot_ = P0/PhysicalConstants::R_J_kmol/T;
    	const double rho_ = cTot_*MW_;
	Product(cTot_, x_, &c_);

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T);
	thermodynamicsMap_.SetPressure(P0);

	// Calculates kinetics
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P0);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_);
	kineticsMap_.FormationRates(&R_);

	for (unsigned int i=1;i<=NC_;++i)	
		S[i] = R_[i]*thermodynamicsMap_.MW()[i];

    	{
		double CpMixMolar; 
		thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
		const double CpMixMass_ = CpMixMolar / MW_;
		const double QR_ = kineticsMap_.HeatRelease(R_);
		
		S[NC_+1]  = QR_;
	}
}

void reactionJacobian( OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap_,
		       const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       Eigen::MatrixXd &J) 
{
     // Calculated as suggested by Buzzi (private communication)
     const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
     const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
     const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
     const double TOLR = 1.e-7;
     const double TOLA = 1.e-12;

     const unsigned int NE = y.Size();

     OpenSMOKE::OpenSMOKEVectorDouble y_plus(NE);
     OpenSMOKE::OpenSMOKEVectorDouble dy_plus(NE);
     OpenSMOKE::OpenSMOKEVectorDouble dy_original(NE);

     for(unsigned int i=1;i<=NE;i++)
		y_plus[i] = y[i];

     // Call equations
     reactionSourceTerms(thermodynamicsMap_, kineticsMap_, y, P0, dy_original);     

     // Derivatives with respect to y[kd]
     for(int kd=1;kd<=NE;kd++)
     {
         double hf = 1.e0;
         double error_weight = 1./(TOLA+TOLR*fabs(y[kd]));
         double hJ = ETA2 * fabs(std::max(y[kd], 1./error_weight));
         double hJf = hf/error_weight;
         hJ = std::max(hJ, hJf);
         hJ = std::max(hJ, ZERO_DER);

         // This is what is done by Buzzi
         double dy = std::min(hJ, 1.e-3 + 1e-3*fabs(y[kd]));
         double udy = 1. / dy;
         y_plus[kd] += dy;

	 reactionSourceTerms(thermodynamicsMap_, kineticsMap_, y_plus, P0, dy_plus);

         for(int j=1;j<=NE;j++)
             J(j-1,kd-1) = (dy_plus[j]-dy_original[j]) * udy;

         y_plus[kd] = y[kd];
     }
 };
