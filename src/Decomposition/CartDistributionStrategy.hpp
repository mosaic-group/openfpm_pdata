#ifndef SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP
#define SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP

#include "Decomposition/AbstractDistributionStrategy.hpp"

template <unsigned int dim, typename domain_type,
          typename AbstractDistStrategy =
              AbstractDistributionStrategy<dim, domain_type>, typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin>
class CartDistributionStrategy {

  using Box = SpaceBox<dim, domain_type>;
  using DGrid = grid_sm<dim, void>;
  using ParmetisGraph = Parmetis<Graph_CSR<nm_v<dim>, nm_e>>;

public:
  CartDistributionStrategy(Vcluster<> &v_cl) : dist(v_cl), graph(v_cl, v_cl.getProcessingUnits()) {
    //! Convert the graph to parmetis => todo rename ParMetisDistribution
  }

  ~CartDistributionStrategy() {}

  void setParameters(DGrid &grid_dec, const Ghost<dim, domain_type> &ghost,
                const grid_sm<dim, void> &sec_dist = grid_sm<dim, void>()) {
    if (sec_dist.size(0) != 0) {
      gr.setDimensions(sec_dist.getSize());
    } else {
      gr = grid_dec;
    }

    dist.ghost = ghost;
  }

  /*! \brief Create the Cartesian graph
   *
   * \param grid info
   * \param dom domain
   */
  void createCartGraph(const size_t (&bc)[dim], ::Box<dim, domain_type> &domain) {
    // Create a cartesian grid graph
    CartesianGraphFactory<dim, Graph_CSR<nm_v<dim>, nm_e>> g_factory_part;
    dist.gp = g_factory_part.template construct<NO_EDGE, nm_v_id, domain_type, dim - 1, 0>(gr.getSize(), domain, bc);
    dist.initLocalToGlobalMap();

    //! Get the number of processing units
    size_t Np = dist.v_cl.getProcessingUnits();

    //! Division of vertices in Np graphs
    //! Put (div+1) vertices in mod graphs
    //! Put div vertices in the rest of the graphs
    size_t mod_v = gr.size() % Np;
    size_t div_v = gr.size() / Np;

    for (size_t i = 0; i <= Np; i++) {
      if (i < mod_v) {
        dist.vtxdist.get(i).id = (div_v + 1) * i;
      } else {
        dist.vtxdist.get(i).id = (div_v)*i + mod_v;
      }
    }

    // Init to 0.0 axis z (todo fix in graphFactory)
    if (dim < 3) {
      for (size_t i = 0; i < dist.gp.getNVertex(); i++) {
        dist.gp.vertex(i).template get<nm_v_x>()[2] = 0.0;
      }
    }

    for (size_t i = 0; i < dist.gp.getNVertex(); i++) {
      dist.gp.vertex(i).template get<nm_v_global_id>() = i;
    }
  }

  void reset() {
    if (dist.is_distributed) {
      graph.reset(dist.gp, dist.vtxdist, dist.m2g, dist.verticesGotWeights);
    } else {
      graph.initSubGraph(dist.gp, dist.vtxdist, dist.m2g, dist.verticesGotWeights);
    }
  }

  void distribute() {
    reset();

    //! Get the processor id
    size_t p_id = dist.v_cl.getProcessUnitID();

    //! Get the number of processing units
    size_t Np = dist.v_cl.getProcessingUnits();

    // Number of local vertex
    size_t nl_vertex = dist.vtxdist.get(p_id + 1).id - dist.vtxdist.get(p_id).id;

    //! Get result partition for this processors
    idx_t *partition = graph.getPartition();

    //! Prepare vector of arrays to contain all partitions
    dist.partitions.get(p_id).resize(nl_vertex);

    if (nl_vertex != 0) {
      std::copy(partition, partition + nl_vertex, &dist.partitions.get(p_id).get(0));
    }

    // Reset data structure to keep trace of new vertices distribution in
    // processors (needed to update main graph)
    for (size_t i = 0; i < Np; ++i) {
      dist.v_per_proc.get(i).clear();
    }

    // Communicate the local distribution to the other processors
    // to reconstruct individually the global graph
    openfpm::vector<size_t> prc;
    openfpm::vector<size_t> sz;
    openfpm::vector<void *> ptr;

    for (size_t i = 0; i < Np; i++) {
      if (i != dist.v_cl.getProcessUnitID()) {
        dist.partitions.get(i).clear();
        prc.add(i);
        sz.add(nl_vertex * sizeof(idx_t));
        ptr.add(dist.partitions.get(p_id).getPointer());
      }
    }

    if (prc.size() == 0) {
      dist.v_cl.sendrecvMultipleMessagesNBX(0, NULL, NULL, NULL, dist.message_receive,
                                       &dist.partitions, NONE);
    } else {
      dist.v_cl.sendrecvMultipleMessagesNBX(prc.size(), &sz.get(0), &prc.get(0),
                                       &ptr.get(0), dist.message_receive,
                                       &dist.partitions, NONE);
    }

    // Update graphs with the received data
    updateGraphs();

    dist.distribute();
  }

  /*! \brief Update main graph ad subgraph with the received data of the
   * partitions from the other processors
   *
   */
  void updateGraphs() {
    dist.sub_sub_owner.clear();

    size_t Np = dist.v_cl.getProcessingUnits();

    // Init n_vtxdist to gather informations about the new decomposition
    openfpm::vector<rid> n_vtxdist(Np + 1);
    for (size_t i = 0; i <= Np; i++)
      n_vtxdist.get(i).id = 0;

    // Update the main graph with received data from processor i
    for (size_t i = 0; i < Np; i++) {
      size_t ndata = dist.partitions.get(i).size();
      size_t k = 0;

      // Update the main graph with the received informations
      for (rid l = dist.vtxdist.get(i); k < ndata && l < dist.vtxdist.get(i + 1);
           k++, ++l) {
        // Create new n_vtxdist (just count processors vertices)
        ++n_vtxdist.get(dist.partitions.get(i).get(k) + 1);

        // vertex id from vtx to grobal id
        auto v_id = dist.m2g.find(l)->second.id;

        // Update proc id in the vertex (using the old map)
        dist.gp.template vertex_p<nm_v_proc_id>(v_id) = dist.partitions.get(i).get(k);

        if (dist.partitions.get(i).get(k) == (long int) dist.v_cl.getProcessUnitID()) {
          dist.sub_sub_owner.add(v_id);
        }

        // Add vertex to temporary structure of distribution (needed to update
        // main graph)
        dist.v_per_proc.get(dist.partitions.get(i).get(k)).add(dist.getVertexGlobalId(l));
      }
    }

    // Create new n_vtxdist (accumulate the counters)
    for (size_t i = 2; i <= Np; i++) {
      n_vtxdist.get(i) += n_vtxdist.get(i - 1);
    }

    // Copy the new decomposition in the main vtxdist
    for (size_t i = 0; i <= Np; i++) {
      dist.vtxdist.get(i) = n_vtxdist.get(i);
    }

    openfpm::vector<size_t> cnt;
    cnt.resize(Np);

    for (size_t i = 0; i < dist.gp.getNVertex(); ++i) {
      size_t pid = dist.gp.template vertex_p<nm_v_proc_id>(i);

      rid j = rid(dist.vtxdist.get(pid).id + cnt.get(pid));
      gid gi = gid(i);

      dist.gp.template vertex_p<nm_v_id>(i) = j.id;
      cnt.get(pid)++;

      dist.setMapId(j, gi);
    }
  }

  ParmetisGraph& getGraph() { return graph; }

  /*! \brief Distribution grid
   *
   * \return the grid
   */
  DGrid getGrid() { return gr; }

// todo private:

  ParmetisGraph graph;

  //! Structure that store the cartesian grid information
  DGrid gr;

  AbstractDistStrategy dist;
};
#endif // SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP
