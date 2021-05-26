//OFstream ofile("output.out", ios_base::app);  // for appending to an existing file
OFstream ofile("output.out");

if(!ofile)      {Info << "File output.out can not be opened\n";return 1;}
