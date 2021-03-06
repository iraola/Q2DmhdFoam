// Steady state termination criteria
epsSteady = readScalar(runTime.controlDict().lookup("epsSteady"));
maxT0 = maxT;
maxT = gMax(T);
repsmaxT = fabs((maxT - maxT0) / (maxT0 + SMALL)) / runTime.deltaT().value();
maxU0 = maxU;
maxU = gMax(mag(U.component(vector::X))());
repsmaxU = fabs((maxU - maxU0) / (maxU0 + SMALL)) / runTime.deltaT().value();

// Write data
if  ((Pstream::master()))
{
    ofile   << runTime.timeName()         << '\t'
            << magUbarStar.value()        << '\t'
            << gradPplus.value()          << '\t'
            << gradP.value()              << '\t'
            << repsmaxU                   << '\t'
            << repsmaxT                 << '\t'
            << endl;
}


scalar maxEpsilon = max(repsmaxT, repsmaxU);
if ((maxEpsilon < epsSteady))
{
    steadyState = true;
    Info << "Steady-state criteria reached" << endl;
    runTime.writeNow();
}
