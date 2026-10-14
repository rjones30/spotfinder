#include <iostream>
#include <TROOT.h>
void makerootvisuals() {
   std::cout << "Building the C++ rootvisuals_C.so and supporting libraries" << std::endl
             << "for the spotfinder tool. Make sure that LD_PRELOAD is not" << std::endl
             << "set in your environment, or you are in for a bad surprise!" << std::endl;
   TString pyinc = gSystem->GetFromPipe(
      "python3 -c \"import sysconfig; print(sysconfig.get_paths()['include'])\"");
   if (pyinc.IsNull() || gSystem->AccessPathName(pyinc + "/pyconfig.h")) {
      std::cerr << "Python headers not found in '" << pyinc
                << "' -- is python3-devel installed?" << std::endl;
      return;
   }
   std::cout << "Using Python headers from " << pyinc << std::endl;
   gSystem->AddIncludePath(("-I" + pyinc).Data());
   gSystem->CompileMacro("Map2D.cc", "kg");
   gSystem->CompileMacro("Couples.C", "kg");
   gSystem->AddIncludePath("-I/cvmfs/oasis.opensciencegrid.org/gluex/Diracxx/x86_64/include");
   gSystem->CompileMacro("CobremsGeneration.cc", "kg");
   gSystem->CompileMacro("rootvisuals.C", "kg");
}
