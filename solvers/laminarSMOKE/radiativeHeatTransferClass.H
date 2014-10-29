
enum radiationModels { RADIATION_NONE, RADIATION_OPTICALLY_THIN} ;

class radiativeHeatTransferClass
{
public:
	radiativeHeatTransferClass()
	{
		radiationModel_ = RADIATION_NONE;
		Tenv_ = 0.;
		Tenv4_ = 0.;

		iCH4_ = 0;
		iCO2_ = 0;
		iCO_  = 0;
		iH2O_ = 0;

		xCH4_ = 0.;
		xCO2_ = 0.;
		xCO_  = 0.;
		xH2O_ = 0.;
	}

	void SetThermodynamics(const OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMapXML)
	{
		iCH4_ = thermodynamicsMapXML.IndexOfSpeciesWithoutError("CH4");
		iH2O_ = thermodynamicsMapXML.IndexOfSpeciesWithoutError("H2O");
		iCO2_ = thermodynamicsMapXML.IndexOfSpeciesWithoutError("CO2");
		iCO_  = thermodynamicsMapXML.IndexOfSpeciesWithoutError("CO");

		if (iCH4_==0) iCH4_ = thermodynamicsMapXML.IndexOfSpeciesWithoutError("ch4");
		if (iH2O_==0) iH2O_ = thermodynamicsMapXML.IndexOfSpeciesWithoutError("h2o");
		if (iCO2_==0) iCO2_ = thermodynamicsMapXML.IndexOfSpeciesWithoutError("co2");
		if (iCO_==0)  iCO_  = thermodynamicsMapXML.IndexOfSpeciesWithoutError("co");
	}

	void SetRadiationModel(const radiationModels model)
	{
		radiationModel_ = model;
	}

	void SetEnvironmentTemperature(const double Tenv)
	{
		Tenv_  = Tenv;
		Tenv4_ = Tenv*Tenv*Tenv*Tenv;
	}

	radiationModels radiationModel() const { return radiationModel_; }

	double aPlanck(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		xCH4_ = (iCH4_>0) ? x[iCH4_] : 0. ;
		xCO2_ = (iCO2_>0) ? x[iCO2_] : 0. ;
		xCO_  = (iCO_>0)  ? x[iCO_]  : 0. ;
		xH2O_ = (iH2O_>0) ? x[iH2O_] : 0. ;

		return aPlanck(T, P_Pa, xCH4_, xH2O_, xCO2_, xCO_);
	}
	
	double Q(const double T, const double a)
	{
		return 4.*sigmaSB_*a*(T*T*T*T - Tenv4_);
	}

private:

	double Tenv_;
	double Tenv4_;
	
	radiationModels radiationModel_;
	
	unsigned int iCH4_;
	unsigned int iCO2_;
	unsigned int iCO_;
	unsigned int iH2O_;

	double xH2O_;
	double xCO_;
	double xCO2_;
	double xCH4_;

	static const double sigmaSB_;	

	double aPlanck(const double T, const double P_Pa, const double xCH4, const double xH2O, const double xCO2, const double xCO);
};

const double radiativeHeatTransferClass::sigmaSB_ = 5.67e-8;		// [W/m2/K4]

double radiativeHeatTransferClass::aPlanck(const double T, const double P_Pa, const double xCH4, const double xH2O, const double xCO2, const double xCO)
{
	double a = 0.;
	const double uT = 1000./T;		// [1]

	// 1. Water [1/m/bar]
	if (iH2O_>0)
	{
		const double KH2O = -0.23093 +uT*(-1.1239+uT*(9.4153 +uT*(-2.9988 +uT*( 0.51382 + uT*(-1.8684e-5)))));
		a+=KH2O*xH2O;
	}

	// 2. Carbon Dioxide [1/m/bar]
	if (iCO2_>0)
	{
		const double KCO2 =  18.741  +uT*(-121.31+uT*(273.5  +uT*(-194.05 +uT*( 56.31   + uT*(-5.8169)))));
		a+=KCO2*xCO2;
	}

	// 3. Carbon monoxide [1/m/bar]
	if (iCO_>0)
	{
		const double KCO = ( T < 750. ) ? 
				4.7869+T*(-0.06953 + T*(2.95775e-4 + T*(-4.25732e-7 + T*2.02894e-10))) :
				10.09+T*(-0.01183 + T*(4.7753e-6 + T*(-5.87209e-10 + T*-2.5334e-14)));
		a+=KCO*xCO;
	}

	// 4. Methane [1/m/bar]
	if (iCH4_>0)
	{
		const double KCH4 = 6.6334 +T*(- 0.0035686+T*(1.6682e-08+T*(2.5611e-10-2.6558e-14*T)));
		a+=KCH4*xCH4;
	}

	// Planck mean absorption coefficient [1/m]
	return a * (P_Pa/1.e5);
}
