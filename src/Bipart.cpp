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


MetisGraph* biparting(PhyDB& db, unsigned Csize, unsigned K) {

  MetisGraph *mG = new MetisGraph();
  GGraph& graph = *(mG->getGraph());
  auto &phy_db_design = *(db.GetDesignPtr());

  int components_count = 0, pins_count = 0, nets_count = 0;
  auto &components = phy_db_design.GetComponentsRef();
  components_count = components.size();
  auto &iopins = phy_db_design.GetIoPinsRef();
  pins_count = iopins.size();
  auto &nets = phy_db_design.GetNetsRef();
  nets_count = nets.size();

  uint32_t hedges = nets_count;
  uint64_t nodes  = components_count;
  //std::cout<<"nodes "<<nodes<<" nets " <<hedges<<"\n";
  if (hedges < 1) { 
    std::cout<<"0 nets\n";
    return mG;
  }

  galois::gstl::Vector<galois::PODResizeableArray<uint32_t>> edges_id(hedges +
                                                                      nodes);
  std::vector<std::vector<EdgeTy>> edges_data(hedges + nodes);
  std::vector<uint64_t> prefix_edges(nodes + hedges);
  std::vector<int> weights(hedges);
  uint32_t edges = 0;
  unsigned nodeid = hedges;
  std::map<std::string, MetisNode> mapNodes;
  std::map<int, MetisNode> mapNodeId;
  uint32_t cnt = 0;
  for (auto &net: nets) {
    std::string net_name(net.GetName());
      //auto &comp_names = net.GetComponentNamesRef();
    //auto &comp_names = net.GetIoPinNamesRef();
    //auto pins = net.GetPinsRef();
    int sz = net.pins_.size();
    if (sz< 1) { 
      std::cout<<"pin name is not defined\n";
    return mG;
  }
    //for (int i = 0; i < comp_names.size(); ++i) {
    //for (auto v : pins) {
    for (auto &pin: net.pins_) {
            std::string component_name =
                phy_db_design.components_[pin.comp_id].GetName();
            std::string macro_name =
                phy_db_design.components_[pin.comp_id].GetMacro()->GetName();
            Macro *macro_ptr = db.GetMacroPtr(macro_name);
            //PhyDBExpects(macro_ptr != nullptr, "Macro does not exist");
            std::string
                pin_name = macro_ptr->GetPinsRef()[pin.pin_id].GetName();
      auto it = mapNodes.find(pin_name);
      if (it != mapNodes.end()) {
        auto valn = it->second;
        int val = valn.nodeid;
        unsigned newval = (val);
        edges_id[cnt].push_back(newval);
        edges++;
      }
      else {
        MetisNode n1;
        n1.name = pin_name;
        n1.nodeid = nodeid++;
        mapNodes[pin_name] = n1;
        mapNodeId[n1.nodeid] = n1;
        edges_id[cnt].push_back(n1.nodeid);
        edges++;
      }
    }
    cnt++;
  }
  uint32_t sizes = hedges + nodes;
  graph.hedges = hedges;
  graph.hnodes = nodes;

  galois::do_all(galois::iterate(uint32_t{0}, (uint32_t) sizes),
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
          graph.getData(n).name = n1.name;
        }
        graph.getData(n).netrand = INT_MAX;
        graph.getData(n).netval  = INT_MAX;
      },
      galois::steal(), galois::loopname("build initial graph"));
  galois::preAlloc(galois::runtime::numPagePoolAllocTotal() * 10);
  galois::reportPageAlloc("MeminfoPre");
  galois::do_all(
      galois::iterate(graph.hedges, graph.size()),
      [&](GNode item) {
        // accum += g->getData(item).getWeight();
        graph.getData(item, galois::MethodFlag::UNPROTECTED)
            .initRefine(0, true);
        graph.getData(item, galois::MethodFlag::UNPROTECTED).initPartition();
      },
      galois::steal(), galois::loopname("initPart"));

  Partitions(mG, Csize, K);

  const int k = K;
  // calculating number of iterations/levels required
  int num = log2(k);

  int kValue[k];
  for (int i = 0; i < k; i++)
    kValue[i] = 0;

  kValue[0]           = (k + 1) / 2;
  kValue[(k + 1) / 2] = k / 2;

  galois::do_all(
      galois::iterate((uint64_t)graph.hedges, (uint64_t)graph.size()),
      [&](GNode n) {
        unsigned pp = graph.getData(n).getPart();
        if (pp == 1) {
          graph.getData(n).setPart((k + 1) / 2);
        }
      },
      galois::steal(), galois::loopname("set part (original graph)"));
  

  // running it level by level

  // toProcess contains nodes to be executed in a given level
  std::set<int> toProcess;
  std::set<int> toProcessNew;
  toProcess.insert(0);
  toProcess.insert((k + 1) / 2);

  std::vector<std::vector<GNode>> nodesvec(k);
  // std::array<std::vector<GNode>, 100> hedgesvec;

  for (int level = 0; level < num; level++) {

    for (int i = 0; i < k; i++)
      nodesvec[i].clear();

    // distributing nodes in relevant vectors according to their current
    // partition assignment
    for (GNode n = graph.hedges; n < graph.size(); n++) {
      unsigned pp = graph.getData(n).getPart();
      nodesvec[pp].push_back(n);
    }

    std::vector<std::vector<GNode>> hedgevec(k);

    // distribute hyperedges according to their current partition
    galois::do_all(
        galois::iterate((uint64_t)0, (uint64_t)graph.hedges),
        [&](GNode h) {
          auto edge = *(graph.edges(h).begin());
          auto dst  = graph.getEdgeDst(edge);
          auto ii   = graph.getData(dst).getPart();

          bool flag = true;

          for (auto n : graph.edges(h)) {
            auto part = graph.getData(graph.getEdgeDst(n)).getPart();

            if (part != ii) {
              flag = false;
              break;
            }
          }

          if (flag)
            graph.getData(h).setPart(ii);
          else
            graph.getData(h).setPart(100000);
        },
        galois::steal(), galois::loopname("distribute hedges"));

    for (GNode h = 0; h < graph.hedges; h++) {
      unsigned part = graph.getData(h).getPart();
      if (part != 100000)
        hedgevec[part].push_back(h);
    }

    // calling Partition for each partition number
    for (unsigned i : toProcess) {
      if (kValue[i] > 1) {
        MetisGraph metisG;
        GGraph& gr = *metisG.getGraph();

        unsigned ed = 0;

        for (auto h : hedgevec[i])
          graph.getData(h).index = ed++;

        unsigned id = ed;
        for (auto n : nodesvec[i]) {
          graph.getData(n).index = id++;
        }

        unsigned totalnodes = id;
        galois::gstl::Vector<galois::PODResizeableArray<uint32_t>> edges_ids(
            totalnodes);
        std::vector<std::vector<EdgeTy>> edge_data(totalnodes);
        std::vector<uint64_t> pre_edges(totalnodes);
        unsigned edges = 0;

        galois::do_all(
            galois::iterate(hedgevec[i]),
            [&](GNode h) {
              for (auto v : graph.edges(h)) {
                auto vv = graph.getEdgeDst(v);

                uint32_t newid = graph.getData(h).index;
                unsigned nm    = graph.getData(vv).index;
                edges_ids[newid].push_back(nm);
              }
            },
            galois::steal(), galois::loopname("populate edge ids"));

        uint64_t num_edges_acc = 0;
        //galois::do_all(
          //  galois::iterate(uint32_t{0}, totalnodes),
            for(uint32_t c = 0;c<totalnodes;c++) {
              pre_edges[c] = edges_ids[c].size();
              num_edges_acc += pre_edges[c];
            }
            //galois::steal(), galois::loopname("set pre edges"));

        edges = num_edges_acc;

        for (uint64_t c = 1; c < totalnodes; ++c) {
          pre_edges[c] += pre_edges[c - 1];
        }
        gr.constructFrom(totalnodes, edges, pre_edges, edges_ids, edge_data);

        gr.hedges = ed;
        gr.hnodes = id - ed;

        galois::do_all(
            galois::iterate(gr),
            [&](GNode n) {
              if (n < gr.hedges)
                gr.getData(n).netnum = n + 1;
              else
                gr.getData(n).netnum = INT_MAX;
              gr.getData(n).netrand = INT_MAX;
              gr.getData(n).netval  = INT_MAX;
              gr.getData(n).nodeid  = n + 1;
            },
            galois::steal(), galois::loopname("build graph: recursion level"));

        Partitions(&metisG, Csize, kValue[i]);

        MetisGraph* mcg = &metisG;

        // now free up the memory by deleting all coarsened graphs
        while (mcg->getCoarserGraph() != NULL) {
          mcg = mcg->getCoarserGraph();
        }

        while (mcg->getFinerGraph() != NULL &&
               mcg->getFinerGraph()->getFinerGraph() != NULL) {
          mcg = mcg->getFinerGraph();
          delete mcg->getCoarserGraph();
        }

        int tmp                   = kValue[i];
        kValue[i]                 = (tmp + 1) / 2;
        kValue[i + (tmp + 1) / 2] = (tmp) / 2;
        toProcessNew.insert(i);
        toProcessNew.insert(i + (tmp + 1) / 2);

        galois::do_all(
            galois::iterate(nodesvec[i]),
            [&](GNode v) {
              GNode n     = graph.getData(v).index;
              unsigned pp = gr.getData(n).getPart();
              if (pp == 0) {
                graph.getData(v).setPart(i);
              } else if (pp == 1) {
                graph.getData(v).setPart(i + (tmp + 1) / 2);
              }
            },
            galois::steal(),
            galois::loopname("set part: inside recursive call"));

        delete mcg;
      } // end if
    }   // end for

    toProcess = toProcessNew;
    toProcessNew.clear();
  } // end while
  return mG;
}
} //namespace
