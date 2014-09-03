{	
	#include "rhoEqn.H"

	while (pimple.loop())
	{
		#include "UEqn.H"
		#include "YEqn.H"
		#include "TEqn.H"

		// Chemical step
		{
			double t0 = runTime.value() - runTime.deltaT().value();
			double tf = runTime.value();
			#include "chemistry.H"
			#include "properties.H"
		}

		// --- Pressure corrector loop
		while (pimple.correct())
		{
			#if OPENFOAM_VERSION == 22
			#include "pEqn.22.H"
			#elif OPENFOAM_VERSION == 23
			#include "pEqn.23.H"
			#endif
		}
	}
}
