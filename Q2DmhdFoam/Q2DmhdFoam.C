/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    Q2DmhdFoam

Description
    Transient solver for 2D buoyancy-driven laminar magnetohydrodynamic
	flow of incompressible conducting fluids using the Boussinesq Model and
	the Sommeria and Moreau (1983) approximation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "readTransportProperties.H"
#   include "readGravitationalAcceleration.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "createTimeControls.H"
#   include "CourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readTimeControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
         ==
           -beta*(T - T0)*g
           -fvm::Sp(mag(B)/tauHa, U)
        );

        solve(UEqn == -fvc::grad(p));
        
        // --- PISO loop

        while (piso.correct())
        {
            volScalarField rUA = 1.0/UEqn.A();

            U = rUA*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf())
                + fvc::ddtPhiCorr(rUA, U, phi);

            adjustPhi(phi, U, p);

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rUA, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve();

                if (piso.finalNonOrthogonalIter())
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            U -= rUA*fvc::grad(p);
            U.correctBoundaryConditions();

            if (magUbar.value()>SMALL)
            {
#               include "pump.H"
            }
        }

        // Solve energy equation
        solve
        (
            fvm::ddt(T)
          + fvm::div(phi, T)
          - fvm::laplacian(DT, T)
          - sourceT/(rho0*Cp)
        );

        rho = rho0*(scalar(1) - beta*(T - T0));

        // Calculate vorticity
        vorticity = fvc::curl(U);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
