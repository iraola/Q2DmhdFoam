// Calculations for average temperature steady state
scalarField m = U.component(vector::X)*mesh.V();
tempVolAverage.value() = gSum(m*T)/gSum(m);
T = T - tempVolAverage;

// Steady state termination criteria
epsSteady = readScalar(runTime.controlDict().lookup("epsSteady"));
maxT0 = maxT;
maxT = gMax(T);
repsmaxT = fabs((maxT - maxT0) / (maxT0 + SMALL)) / runTime.deltaT().value();
maxU0 = maxU;
maxU = gMax(mag(U.component(vector::X))());
repsmaxU = fabs((maxU - maxU0) / (maxU0 + SMALL)) / runTime.deltaT().value();

scalar maxEpsilon = max(repsmaxT, repsmaxU);
if ((maxEpsilon < epsSteady))
{
    steadyState = true;
    Info << "Steady-state criteria reached" << endl;
    runTime.writeNow();
}
