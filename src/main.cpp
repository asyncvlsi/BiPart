#include <iostream>
#include "Bipart.h"
//using namespace bipart;

using namespace phydb;
int main(int argc, char** argv) {

  PhyDB db;

  string lefFileName = "benchmark_1K.lef";
  db.ReadLef(lefFileName);

  string defFileName = "benchmark_1K.def";
  db.ReadDef(defFileName);

  //auto &phy_db_design = *(db.GetDesignPtr());
  int Csize = 25;
  int Rsize = 2;
  bipart::biparting(db, Csize, Rsize);
  return 0;
};
