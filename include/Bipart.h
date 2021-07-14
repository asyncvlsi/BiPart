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

//#include "/net/ohm/export/iss/sepideh/phyDB/include/phydb.h"
#include "phydb/phydb.h"

namespace bipart {
void Partition(MetisGraph* metisGraph, unsigned coarsenTo, unsigned K);

int hash(unsigned val);
MetisGraph* biparting(phydb::PhyDB& db, unsigned Csize, unsigned K);

}
