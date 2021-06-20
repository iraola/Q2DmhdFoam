// Calculations for average temperature steady state
scalarField m = U.component(vector::X)*mesh.V();
tempVolAverage.value() = gSum(m*T)/gSum(m);
T = T - tempVolAverage;
