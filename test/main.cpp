#include <iostream>
#include "Bipart.h"
#include "Metis.h"
#include <string>
//#include "phydb/phydb.h"
//#include <phydb/phydb.h>
//#include "/net/ohm/export/iss/sepideh/phyDB/include/phydb.h"
static const char* name = "BIPART";
static const char* desc =
    "Partitions a hypergraph into K parts and minimizing the graph cut";
static const char* url = "BiPart";
using namespace phydb;
using namespace bipart;
//using namespace std;
int main(int argc, char** argv) {
  galois::SharedMemSys G;
  //std::string inputFile = "benchmark_1K.lef";
  //LonestarStart(argc, argv, name, desc, url, &inputFile);


  PhyDB db;

  std::string lefFileName = argv[1];//"benchmark_1K.lef";
  db.ReadLef(lefFileName);

  std::string defFileName = argv[2];//"benchmark_1K.def";
  db.ReadDef(defFileName);

  int Csize = std::stoi(argv[3]);//Coarsening levels;
  int K = std::stoi(argv[4]);//number of partitions;
  std::cout << "reading lef/def done\n";
  MetisGraph* metisgraph = biparting(db, Csize, K);
  return 0;
};
