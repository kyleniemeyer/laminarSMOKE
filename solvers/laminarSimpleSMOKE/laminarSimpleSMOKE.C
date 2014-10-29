/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laminarSMOKE

Description
    Solver for reactive flows with detailed kinetic mechanisms.

\*---------------------------------------------------------------------------*/

// This is a steady state simulation
#define STEADYSTATE 1

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"

// Reactor utilities
#include "reactors/utilities/Utilities"

// OpenFOAM
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"

// Additional include files
#include "sparkModel.H"
#include "radiativeHeatTransferClass.H"
#include "utilities.H"

// Linearization
#include "linearizationFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "readGravitationalAcceleration.H"

	simpleControl simple(mesh);

	#include "createBasicFields.H"
	#include "readOptions.H"
	#include "createChemicalFields.H"
	#include "createFvOptions.H"
	#include "memoryAllocation.H"
	#include "properties.H"
	#include "createAdditionalFields.H"
	#include "initContinuityErrs.H"

	dimensionedScalar initialMass = fvc::domainIntegrate(rho);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
         Info<< "Time = " << runTime.timeName() << nl << endl;
		
	 double t0 = runTime.value() - runTime.deltaT().value();
	 double tf = runTime.value();

	 // Pressure-velocity SIMPLE corrector
         {
            #include "UEqn.H"

	    #include "properties.H"
	    #include "jacobianEvaluation.H"
            #include "YEqn.H"
	    #include "TEqn.H"
            #include "pEqn.H"
         }	
				
	 runTime.write();
		
         Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
              << "  ClockTime = " << runTime.elapsedClockTime() << " s"
              << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
