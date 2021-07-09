#include "Metis.h"
#include <galois/graphs/ReadGraph.h>
#include <galois/Timer.h>
//#include <Lonestar/BoilerPlate.h>
#include <galois/graphs/FileGraph.h>
#include <galois/LargeArray.h>
#include <galois/gstl.h>

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <array>
#include <unordered_set>

#include "/net/ohm/export/iss/sepideh/phyDB/include/phydb.h"

//namespace cll = llvm::cl;
/*extern cll::opt<float> tolerance;
                                 
extern cll::opt<std::string>inputFile;

extern cll::opt<bool>verbose;
extern cll::opt<std::string>outfile;
extern cll::opt<unsigned> csize;
extern cll::opt<unsigned> refiter;
extern cll::opt<unsigned> numPartitions;
extern cll::opt<scheduleMode> schedulingMode;
extern cll::opt<bool>hmetisgraph;
extern cll::opt<bool>output;
extern cll::opt<double> imbalance;*/

void Partition(MetisGraph* metisGraph, unsigned coarsenTo, unsigned K);

int hash(unsigned val);
void biparting(int argc, char** argv);
