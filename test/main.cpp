#include <iostream>
#include "Bipart.h"
#include "Metis.h"

using namespace phydb;
using namespace bipart;
int main(int argc, char** argv) {

  PhyDB db;

  string lefFileName = argv[1];//"benchmark_1K.lef";
  db.ReadLef(lefFileName);

  string defFileName = argv[2];//"benchmark_1K.def";
  db.ReadDef(defFileName);

  int Csize = stoi(argv[3]);//Coarsening levels;
  int K = stoi(argv[4]);//number of partitions;

  MetisGraph* metisgraph = biparting(db, Csize, K);
  return 0;
};
