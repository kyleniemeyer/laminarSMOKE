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
			#include "pEqn.H"
		}
	}
}
