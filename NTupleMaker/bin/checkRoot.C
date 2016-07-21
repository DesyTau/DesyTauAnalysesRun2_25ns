#include "TFile.h"

int main(int argc, const char** argv){

  if (argc != 2)
    return 1;
  
  TFile* f = TFile::Open(argv[1], "read");

  if (f==0)
    return 2;

  if (f->GetListOfKeys()->GetSize() == 0)
    return 3; 

  if (f->GetEND() > f->GetSize())
    return 4; 

  if (f->GetSeekKeys()<=f->GetEND()-f->GetSize())
    return 5;
  
  return 0;
}
