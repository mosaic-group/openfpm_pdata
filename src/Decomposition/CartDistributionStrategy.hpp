#ifndef SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP
#define SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP

#include "Decomposition/Domain_NN_calculator_cart.hpp"
#include "Decomposition/AbstractDistributionStrategy.hpp"

template <unsigned int dim, typename domain_type, typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin, DecompositionGraph = Graph_CSR<nm_v<dim>, nm_e>,
          typename AbstractDistStrategy =
              AbstractDistributionStrategy<dim, domain_type, Parmetis<DecompositionGraph>>>
class CartDistributionStrategy : public domain_nn_calculator_cart<dim> {

public:
  CartDistributionStrategy(Vcluster<> &v_cl) : dist(v_cl), graph(v_cl, v_cl.getProcessingUnits()) {
    //! Convert the graph to parmetis => todo rename ParMetisDistribution
  }

  ~CartDistributionStrategy() {}

  void setParameters(grid_sm<dim, void> &gr, const Ghost<dim, domain_type> &ghost, const grid_sm<dim, void> &sec_dist = grid_sm<dim, void>(), const size_t (&bc)[dim], ::Box<dim, domain_type> &box) {
    if (sec_dist.size(0) != 0) {
      gr_dist.setDimensions(sec_dist.getSize());
    } else {
      gr_dist = gr;
    }

    inner.setGhost(ghost);
    createCartGraph(bc, box);
  }

  /*! \brief Create the Cartesian graph
   *
   * \param grid info
   * \param dom domain
   */
  void createCartGraph(DecompositionGraph& gp, const size_t (&bc)[dim], ::Box<dim, domain_type> &domain) {
    // Create a cartesian grid graph
    CartesianGraphFactory<dim, DecompositionGraph> g_factory_part;
    gp = g_factory_part.template construct<NO_EDGE, nm_v_id, domain_type, dim - 1, 0>(gr_dist.getSize(), domain, bc);
    inner.initLocalToGlobalMap();

    //! Get the number of processing units
    size_t Np = inner.getVcluster().getProcessingUnits();

    //! Division of vertices in Np graphs
    //! Put (div+1) vertices in mod graphs
    //! Put div vertices in the rest of the graphs
    size_t mod_v = gr_dist.size() % Np;
    size_t div_v = gr_dist.size() / Np;

    for (size_t i = 0; i <= Np; i++) {
      if (i < mod_v) {
        inner.getVtxdist().get(i).id = (div_v + 1) * i;
      } else {
        inner.getVtxdist().get(i).id = (div_v)*i + mod_v;
      }
    }

    // Init to 0.0 axis z (todo fix in graphFactory)
    if (dim < 3) {
      for (size_t i = 0; i < gp.getNVertex(); i++) {
        gp.vertex(i).template get<nm_v_x>()[2] = 0.0;
      }
    }

    for (size_t i = 0; i < gp.getNVertex(); i++) {
      gp.vertex(i).template get<nm_v_global_id>() = i;
    }
  }

  /*! \brief Return the global id of the owned sub-sub-domain
   *
   * \param id in the list of owned sub-sub-domains
   *
   * \return the global id
   *
   */
  size_t getOwnerSubSubDomain(size_t id) const { return sub_sub_owner.get(id); }

  /*! \brief Return the total number of sub-sub-domains this processor own
   *
   * \return the total number of sub-sub-domains owned by this processor
   *
   */
  size_t getNOwnerSubSubDomains() const { return sub_sub_owner.size(); }

  /*! \brief function that return the position of the vertex in the space
   *
   * \param id vertex id
   * \param pos vector that will contain x, y, z
   *
   */
  void getSubSubDomainPosition(size_t id, T (&pos)[dim]) {
    // Copy the geometrical informations inside the pos vector
    pos[0] = gp.vertex(id).template get<nm_v_x>()[0];
    pos[1] = gp.vertex(id).template get<nm_v_x>()[1];
    if (dim == 3) {
      pos[2] = gp.vertex(id).template get<nm_v_x>()[2];
    }
  }

  void distribute(DecompositionGraph& gp) {
    reset(gp);
    graph.decompose(inner.getVtxdist());  // this is the distributing step
    postDecomposition();  // todo rename
  }

  void postDecomposition() {
    //! Get the processor id
    size_t p_id = inner.getVcluster().getProcessUnitID();

    //! Get the number of processing units
    size_t Np = inner.getVcluster().getProcessingUnits();

    // Number of local vertex
    size_t nl_vertex = inner.getVtxdist().get(p_id + 1).id - inner.getVtxdist().get(p_id).id;

    //! Get result partition for this processors
    idx_t *partition = graph.getPartition();

    //! Prepare vector of arrays to contain all partitions
    inner.partitions.get(p_id).resize(nl_vertex);

    if (nl_vertex != 0) {
      std::copy(partition, partition + nl_vertex, &inner.partitions.get(p_id).get(0));
    }

    // Reset data structure to keep trace of new vertices distribution in
    // processors (needed to update main graph)
    for (size_t i = 0; i < Np; ++i) {
      inner.v_per_proc.get(i).clear();
    }

    // Communicate the local distribution to the other processors
    // to reconstruct individually the global graph
    openfpm::vector<size_t> prc;
    openfpm::vector<size_t> sz;
    openfpm::vector<void *> ptr;

    for (size_t i = 0; i < Np; i++) {
      if (i != inner.getVcluster().getProcessUnitID()) {
        inner.partitions.get(i).clear();
        prc.add(i);
        sz.add(nl_vertex * sizeof(idx_t));
        ptr.add(inner.partitions.get(p_id).getPointer());
      }
    }

    if (prc.size() == 0) {
      inner.getVcluster().sendrecvMultipleMessagesNBX(0, NULL, NULL, NULL, inner.message_receive,
                                       &inner.partitions, NONE);
    } else {
      inner.getVcluster().sendrecvMultipleMessagesNBX(prc.size(), &sz.get(0), &prc.get(0),
                                       &ptr.get(0), inner.message_receive,
                                       &inner.partitions, NONE);
    }

    // Update graphs with the received data
    updateGraphs();

    inner.distribute();
  }

  /*! \brief Refine current decomposition
   *
   * It makes a refinement of the current decomposition using Parmetis function
   * RefineKWay After that it also does the remapping of the graph
   */
  void refine(DecompositionGraph& gp) {
    graph.reset(gp, vtxdist, m2g, verticesGotWeights); // reset
    graph.refine(inner.getVtxdist());                             // refine
    distribute();
  }

  void onEnd() {
    domain_nn_calculator_cart<dim>::reset();
    domain_nn_calculator_cart<dim>::setParameters(proc_box);
  }

  void reset(DecompositionGraph& gp) {
    if (inner.is_distributed) {
      graph.reset(gp, inner.getVtxdist(), inner.m2g, inner.verticesGotWeights);
    } else {
      graph.initSubGraph(gp, inner.getVtxdist(), inner.m2g, inner.verticesGotWeights);
    }
  }

  /*! \brief operator to init ids vector
   *
   * operator to init ids vector
   *
   */
  void initLocalToGlobalMap() {
    gid g;
    rid i;
    i.id = 0;

    m2g.clear();
    for (; (size_t)i.id < gp.getNVertex(); ++i) {
      g.id = i.id;

      m2g.insert({i, g});
    }
  }

  /*! \brief Update main graph ad subgraph with the received data of the
   * partitions from the other processors
   *
   */
  void updateGraphs() {
    inner.sub_sub_owner.clear();

    size_t Np = inner.getVcluster().getProcessingUnits();

    // Init n_vtxdist to gather informations about the new decomposition
    openfpm::vector<rid> n_vtxdist(Np + 1);
    for (size_t i = 0; i <= Np; i++)
      n_vtxdist.get(i).id = 0;

    // Update the main graph with received data from processor i
    for (size_t i = 0; i < Np; i++) {
      size_t ndata = inner.partitions.get(i).size();
      size_t k = 0;

      // Update the main graph with the received informations
      for (rid l = inner.getVtxdist().get(i); k < ndata && l < inner.getVtxdist().get(i + 1);
           k++, ++l) {
        // Create new n_vtxdist (just count processors vertices)
        ++n_vtxdist.get(inner.partitions.get(i).get(k) + 1);

        // vertex id from vtx to grobal id
        auto v_id = inner.m2g.find(l)->second.id;

        // Update proc id in the vertex (using the old map)
        gp.template vertex_p<nm_v_proc_id>(v_id) = inner.partitions.get(i).get(k);

        if (inner.partitions.get(i).get(k) == (long int) inner.getVcluster().getProcessUnitID()) {
          inner.sub_sub_owner.add(v_id);
        }

        // Add vertex to temporary structure of distribution (needed to update
        // main graph)
        inner.v_per_proc.get(inner.partitions.get(i).get(k)).add(inner.getVertexGlobalId(l));
      }
    }

    // Create new n_vtxdist (accumulate the counters)
    for (size_t i = 2; i <= Np; i++) {
      n_vtxdist.get(i) += n_vtxdist.get(i - 1);
    }

    // Copy the new decomposition in the main vtxdist
    for (size_t i = 0; i <= Np; i++) {
      inner.getVtxdist().get(i) = n_vtxdist.get(i);
    }

    openfpm::vector<size_t> cnt;
    cnt.resize(Np);

    for (size_t i = 0; i < gp.getNVertex(); ++i) {
      size_t pid = gp.template vertex_p<nm_v_proc_id>(i);

      rid j = rid(inner.getVtxdist().get(pid).id + cnt.get(pid));
      gid gi = gid(i);

      gp.template vertex_p<nm_v_id>(i) = j.id;
      cnt.get(pid)++;

      inner.setMapId(j, gi);
    }
  }

private:

  //! Id of the sub-sub-domain where we set the costs
  openfpm::vector<size_t> sub_sub_owner;

  //! Processor domain bounding box
  ::Box<dim, size_t> proc_box;

  //! Structure that store the cartesian grid information
  grid_sm<dim, void> gr;

  AbstractDistStrategy inner;
};
#endif // SRC_DECOMPOSITION_CART_DISTRIBUTION_STRATEGY_HPP
