#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

int main (int argc, char * argv[]) {

  using namespace std;

  if (argc!=3) {
    std::cout << "Usage : Splitter <SampleName> <nperiod>" << std::endl;
    exit(-1);
  }

  std::ifstream input0(argv[1]);
  std::ifstream input(argv[1]);
  std::string sampleName(argv[1]);

  int period = atoi(argv[2]);

  int nfiles = 0;
  std::string dummy;
  while (input0>>dummy) nfiles++;
  
  unsigned int ncycles = nfiles / period;
  unsigned int nreminder = nfiles % period;

  std::cout << std::endl;
  std::cout << "number of files  = " << nfiles << std::endl;
  std::cout << "period           = " << period << std::endl;
  std::cout << "number of cycles = " << ncycles << std::endl;
  std::cout << "reminder         = " << nreminder << std::endl;
  std::cout << std::endl;

  for (unsigned int i=0; i<ncycles; ++i) {
    char outputNumber[10];
    if (i<10)
      sprintf(outputNumber,"%1i",i);
    else if (i<100)
      sprintf(outputNumber,"%2i",i);
    else 
      sprintf(outputNumber,"%3i",i);
    std::string outputName = sampleName + "_files/" + sampleName + "_" + std::string(outputNumber);
    std::ofstream outputFile(outputName);
    for (int j=0; j<period; ++j) {
      input >> dummy;
      //      std::cout << dummy << std::endl;
      outputFile << dummy << std::endl;
    }

  }
  char outputNumber[10];
  if (ncycles<10)
    sprintf(outputNumber,"%1i",ncycles);
  else if (ncycles<100)
    sprintf(outputNumber,"%2i",ncycles);
  else
    sprintf(outputNumber,"%3i",ncycles);
  std::string outputName = sampleName + "_files/" + sampleName + "_" + std::string(outputNumber);
  std::ofstream outputFile(outputName);
  for (unsigned int j=0; j<nreminder; ++j) {
    input >> dummy;
    //    std::cout << dummy << std::endl;
    outputFile << dummy << std::endl;
  }


}
