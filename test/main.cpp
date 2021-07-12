#include <iostream>
#include "Bipart.h"

using namespace phydb;
int main(int argc, char** argv) {

  PhyDB db;

  string lefFileName = argv[1];//"benchmark_1K.lef";
  db.ReadLef(lefFileName);

  string defFileName = argv[2];//"benchmark_1K.def";
  db.ReadDef(defFileName);

  int Csize = stoi(argv[3]);//25;
  int Rsize = stoi(argv[4]);//2;

  bipart::biparting(db, Csize, Rsize);
  return 0;
};
