#include "Bipart.h"
using namespace phydb;
//namespace cll = llvm::cl;

// cll::opt<float> tolerance("tolerance", cll::desc("tolerance"),
 //                                cll::init(tolerance));
// const char* name = "bipart";
 //const char* desc =
   // "partitions a hypergraph into k parts and minimizing the graph cut";
 //const char* url = "hyperbip";

 /*cll::opt<std::string>
    inputFile(cll::Positional, cll::desc("<input file>"));

 cll::opt<bool> weighted("weighted", cll::desc("weighted"),
                               cll::init(false));
 cll::opt<bool>
    verbose("verbose",
            cll::desc("verbose output (debugging mode, takes extra time)"),
            cll::init(false));
 cll::opt<std::string> outfile("outputfile",
                                     cll::desc("output partition file name"));
 cll::opt<unsigned> csize(cll::Positional,
                                cll::desc("<size of coarsest graph>"),
                                cll::init(25));

 cll::opt<unsigned> refiter(cll::Positional,
                                  cll::desc("<number of iterations in ref>"),
                                  cll::init(2));
 cll::opt<unsigned> numPartitions(cll::Positional,
                                        cll::desc("<number of partitions>"),
                                        cll::init(2));

 cll::opt<scheduleMode> schedulingmode(
    cll::desc("choose a inital scheduling mode:"),
    cll::values(clEnumVal(PLD, "pld"), clEnumVal(PP, "pp"), clEnumVal(WD, "wd"),
                clEnumVal(RI, "ri"), clEnumVal(MRI, "mri"),
                clEnumVal(MDEG, "mdeg"), clEnumVal(DEG, "deg"),
                clEnumVal(MWD, "mwd"), clEnumVal(HIS, "his"),
                clEnumVal(RAND, "random")),
    cll::init(RAND));


//! flag that forces user to be aware that they should be passing in a
//! hmetis graph.
 cll::opt<bool>
    hmetisgraph("hmetisgraph",
                cll::desc("specify that the input graph is a hmetis"),
                cll::init(false));

 cll::opt<bool>
    output("output", cll::desc("specify if partitions need to be written"),
           cll::init(false));


 cll::opt<double> imbalance(
    "balance",
    cll::desc("percentage deviated from mean partition size (default 5)"),
    cll::init(5.0));*/
double Ptime = 0.0f;
double Ctime = 0.0f;
double Rtime = 0.0f;

constexpr  const unsigned chunk_size = 32;


/**
 * partitioning
 */

void Partition(MetisGraph* metisGraph, unsigned coarsenTo, unsigned K) {
  galois::StatTimer execTime("Timer_0");
  execTime.start();

  galois::StatTimer T("CoarsenSEP");
  scheduleMode schedulingMode = RAND;
  T.start();
  MetisGraph* mcg = coarsen(metisGraph, coarsenTo, schedulingMode);
  T.stop();
  std::cout<<"end of coarsre\n";

  galois::StatTimer T2("PartitionSEP");
  T2.start();
  partition(mcg, K);
  std::cout<<"end of part\n";
  T2.stop();

  galois::StatTimer T3("Refine");
  T3.start();
  refine(mcg, K, 5.0);
  T3.stop();
  std::cout<<"end of refine\n";
  Ctime += (T.get()/1000.0f);
  Ptime += (T2.get()/1000.0f);
  Rtime += (T3.get()/1000.0f);

  execTime.stop();
}



int hash(unsigned val) {
  unsigned long int seed = val * 1103515245 + 12345;
  return ((unsigned)(seed / 65536) % 32768);
}

void biparting(int argc, char** argv) {
  galois::SharedMemSys G;
  //std::string inputFile = "input.txt".c_str();
  //LonestarStart(argc, argv, name, desc, url, &inputFile);

  galois::StatTimer totalTime("TimerTotal");
  totalTime.start();

  MetisGraph metisGraph;
  GGraph& graph = *metisGraph.getGraph();


  PhyDB db;

    string lefFileName = "benchmark_1K.lef";
    db.ReadLef(lefFileName);

    string defFileName = "benchmark_1K.def";
    db.ReadDef(defFileName);
    std::cout<<"after reading lef and def\n";

    auto &phy_db_design = *(db.GetDesignPtr());

    int components_count = 0, pins_count = 0, nets_count = 0;
    auto &components = phy_db_design.GetComponentsRef();
    components_count = components.size();
    auto &iopins = phy_db_design.GetIoPinsRef();
    pins_count = iopins.size();
    auto &nets = phy_db_design.GetNetsRef();
    nets_count = nets.size();
    std::cout<<nets_count<<" "<<pins_count<<" "<<components_count<<"\n";
    auto die_area = phy_db_design.GetDieArea();
    int LLX = die_area.LLX();
    int LLY = die_area.LLX();
    std::cout<<LLX<<"\n";
    std::cout<<LLY<<"\n";


  uint32_t hedges = nets_count;
  uint64_t nodes  = components_count;
  std::cout << "hedges: " << hedges << "\n";
  std::cout << "nodes: " << nodes << "\n\n";

  galois::StatTimer T("buildingG");
  T.start();

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
      //auto &pin_names = net.GetPinNamesRef();
      //auto &iopin_names = net.GetIoPinNamesRef();
      //int net_capacity = comp_names.size();
      int sz = comp_names.size();
      for (int i = 0; i < sz; ++i) {
        auto valn = mapNodes[comp_names[i]];
        int val = valn.nodeid;
        unsigned newval = (val);
        edges_id[cnt].push_back(newval);
        edges++;
        cnt++;
        //AddBlkPinToNet(comp_names[i], pin_names[i], net_name);
      }
    }
  uint32_t sizes = hedges + nodes;
  graph.hedges = hedges;
  graph.hnodes = nodes;
  std::cout<<"edges are "<<edges<<"\n";

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
          //graph.getData(n).setWeight(weights[n]); 
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
  T.stop();
  std::cout << "time to build a graph " << T.get() << "\n";
  graphStat(graph);
  std::cout << "\n";
  galois::preAlloc(galois::runtime::numPagePoolAllocTotal() * 10);
  galois::reportPageAlloc("MeminfoPre");
Partition(&metisGraph, 25, 2);

  
std::cout<<edges<<"\n";
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
