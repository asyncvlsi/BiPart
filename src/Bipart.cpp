#include "Bipart.h"
using namespace phydb;
namespace bipart {

double Ptime = 0.0f;
double Ctime = 0.0f;
double Rtime = 0.0f;

constexpr const unsigned chunk_size = 32;


/**
 * partitioning
 */

void Partitions(MetisGraph* metisGraph, unsigned coarsenTo, unsigned K) {
  galois::StatTimer execTime("Timer_0");
  execTime.start();

  galois::StatTimer T("CoarsenSEP");
  scheduleMode schedulingMode = RAND;
  T.start();
  MetisGraph* mcg = coarsen(metisGraph, coarsenTo, schedulingMode);
  T.stop();

  galois::StatTimer T2("PartitionSEP");
  T2.start();
  partition(mcg, K);
  T2.stop();

  galois::StatTimer T3("Refine");
  T3.start();
  refine(mcg, K, 5.0);
  T3.stop();
  Ctime += (T.get()/1000.0f);
  Ptime += (T2.get()/1000.0f);
  Rtime += (T3.get()/1000.0f);

  execTime.stop();
}


void biparting(PhyDB& db, int Csize, int Rsize) {
  galois::SharedMemSys G;


  MetisGraph metisGraph;
  GGraph& graph = *metisGraph.getGraph();

  auto &phy_db_design = *(db.GetDesignPtr());

  int components_count = 0, pins_count = 0, nets_count = 0;
  auto &components = phy_db_design.GetComponentsRef();
  components_count = components.size();
  auto &iopins = phy_db_design.GetIoPinsRef();
  pins_count = iopins.size();
  auto &nets = phy_db_design.GetNetsRef();
  nets_count = nets.size();
  auto die_area = phy_db_design.GetDieArea();
  int LLX = die_area.LLX();
  int LLY = die_area.LLX();

  uint32_t hedges = nets_count;
  uint64_t nodes  = components_count;
  std::cout << "Nets: " << hedges << "\n";
  std::cout << "std cells: " << nodes << "\n\n";

  galois::gstl::Vector<galois::PODResizeableArray<uint32_t>> edges_id(hedges +
                                                                      nodes);
  std::vector<std::vector<EdgeTy>> edges_data(hedges + nodes);
  std::vector<uint64_t> prefix_edges(nodes + hedges);
  std::vector<int> weights(hedges);
  uint32_t edges = 0;
    unsigned nodeid = hedges;
    std::map<std::string, MetisNode> mapNodes;
    std::map<int, MetisNode> mapNodeId;
    for (auto &comp: components) {
      std::string blk_name(comp.GetName());
      std::string blk_type_name(comp.GetMacroName());
      auto location = comp.GetLocation();
      int llx = location.x;
      int lly = location.y;
     // GNode node;
      MetisNode n1;
      n1.name = blk_name;
      n1.area = llx*lly;
      n1.nodeid = nodeid++;
      mapNodes[blk_name] = n1;
      mapNodeId[n1.nodeid] = n1;
    }
    uint32_t cnt = 0;
    for (auto &net: nets) {
      std::string net_name(net.GetName());
      auto &comp_names = net.GetComponentNamesRef();
      int sz = comp_names.size();
      for (int i = 0; i < sz; ++i) {
        auto valn = mapNodes[comp_names[i]];
        int val = valn.nodeid;
        unsigned newval = (val);
        edges_id[cnt].push_back(newval);
        edges++;
        cnt++;
      }
    }
  uint32_t sizes = hedges + nodes;
  graph.hedges = hedges;
  graph.hnodes = nodes;

  galois::do_all(galois::iterate(uint32_t{0}, sizes),
                 [&](uint32_t c) { prefix_edges[c] = edges_id[c].size(); });

  for (uint64_t c = 1; c < nodes + hedges; ++c) {
    prefix_edges[c] += prefix_edges[c - 1];
  }

  graph.constructFrom(nodes + hedges, edges, prefix_edges, edges_id,
                      edges_data);
  galois::do_all(
      galois::iterate(graph),
      [&](GNode n) {
        if (n < hedges) {
          graph.getData(n).netnum = n + 1;
          graph.getData(n).nodeid  = n + 1;
        }
        else {
          graph.getData(n).netnum = INT_MAX;
          MetisNode n1 = mapNodeId[n];
          graph.getData(n).nodeid  = n + 1;
          graph.getData(n).area = n1.area;
          graph.getData(n).name = n1.name;
        }
        graph.getData(n).netrand = INT_MAX;
        graph.getData(n).netval  = INT_MAX;
      },
      galois::steal(), galois::loopname("build initial graph"));
  galois::preAlloc(galois::runtime::numPagePoolAllocTotal() * 10);
  galois::reportPageAlloc("MeminfoPre");
  Partitions(&metisGraph, Csize, Rsize);

  
  galois::GAccumulator<unsigned int> area0;
  galois::GAccumulator<unsigned int> area1;
  
  galois::do_all(
      galois::iterate(graph),
      [&](GNode n) {
        if (n > hedges) {
          int part = graph.getData(n).getPart();
          if (part == 0) 
            area0 += graph.getData(n).area;
          else
            area1 += graph.getData(n).area;
        }
      },
      galois::steal(), galois::loopname("build initial graph"));

  int part1x = LLX/2;
  int part1y = LLY/2;
}
} //namespace
