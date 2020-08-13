#ifndef SRC_DECOMPOSITION_ABSTRACTDISTRIBUTIONSTRATEGY_HPP
#define SRC_DECOMPOSITION_ABSTRACTDISTRIBUTIONSTRATEGY_HPP

#include "Space/Ghost.hpp"
#include "Decomposition/Domain_NN_calculator_cart.hpp"
#include "Graph/CartesianGraphFactory.hpp"
#include "Graph/ids.hpp"
#include "SubdomainGraphNodes.hpp"

/*! \brief Class that distribute sub-sub-domains across processors
 */
template <unsigned int dim,
          typename T,
          typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin,
          typename DGrid = grid_sm<dim, void>>
class AbstractDistributionStrategy : public domain_nn_calculator_cart<dim> {
  //! It simplify to access the SpaceBox element
  using Box = SpaceBox<dim, T>;

public:
  //! Vcluster
  Vcluster<>& v_cl;

  /*! Constructor
   *
   * \param v_cl Vcluster to use as communication object in this class
   */
  AbstractDistributionStrategy(Vcluster<>& v_cl)
    : is_distributed(false),
      v_cl(v_cl),
      vtxdist(v_cl.getProcessingUnits() + 1),
      partitions(v_cl.getProcessingUnits()),
      v_per_proc(v_cl.getProcessingUnits()) {}

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

  /*! \brief Returns total number of sub-sub-domains in the distribution graph
   *
   * \return the total number of sub-sub-domains
   *
   */
  size_t getNSubSubDomains() const { return gp.getNVertex(); }

  /*! \brief function that return the position of the vertex in the space
   *
   * \param id vertex id
   * \param pos vector that will contain x, y, z
   *
   */
  void getSubSubDomainPosition(size_t id, T (&pos)[dim]) {
#ifdef SE_CLASS1
    if (id >= gp.getNVertex())
      std::cerr << __FILE__ << ":" << __LINE__
                << "Such vertex doesn't exist (id = " << id << ", "
                << "total size = " << gp.getNVertex() << ")\n";
#endif

    // Copy the geometrical informations inside the pos vector
    pos[0] = gp.vertex(id).template get<nm_v_x>()[0];
    pos[1] = gp.vertex(id).template get<nm_v_x>()[1];
    if (dim == 3) {
      pos[2] = gp.vertex(id).template get<nm_v_x>()[2];
    }
  }

  /*! \brief Set migration cost of the vertex id
   *
   * \param id of the vertex to update
   * \param migration cost of the migration
   */
  void setMigrationCost(size_t id, size_t migration) {
#ifdef SE_CLASS1
    if (id >= gp.getNVertex())
      std::cerr << __FILE__ << ":" << __LINE__
                << "Such vertex doesn't exist (id = " << id << ", "
                << "total size = " << gp.getNVertex() << ")\n";
#endif

    gp.vertex(id).template get<nm_v_migration>() = migration;
  }

  /*! \brief function that get the weight of the vertex
   * (computation cost of the sub-sub-domain id)
   *
   * \param id vertex id
   *
   */
  size_t getSubSubDomainComputationCost(size_t id) {
#ifdef SE_CLASS1
    if (id >= gp.getNVertex())
      std::cerr << __FILE__ << ":" << __LINE__
                << "Such vertex doesn't exist (id = " << id << ", "
                << "total size = " << gp.getNVertex() << ")\n";
#endif

    return gp.vertex(id).template get<nm_v_computation>();
  }

  /*! \brief Add computation cost i to the subsubdomain with global id gid
   *
   * \param gid global id of the subsubdomain to update
   * \param i Cost increment
   */
  void addComputationCost(size_t gid, size_t i) {
    size_t c = getSubSubDomainComputationCost(gid);
    setComputationCost(gid, c + i);
  }

  /*! \brief Set communication cost of the edge id
   *
   * \param v_id Id of the source vertex of the edge
   * \param e i child of the vertex
   * \param communication Communication value
   */
  void setCommunicationCost(size_t v_id, size_t e, size_t communication) {
#ifdef SE_CLASS1

    size_t e_id = v_id + e;

    if (e_id >= gp.getNEdge())
      std::cerr << "Such edge doesn't exist (id = " << e_id << ", "
                << "total size = " << gp.getNEdge() << ")\n";
#endif

    gp.getChildEdge(v_id, e).template get<nm_e::communication>() =
        communication;
  }

  /*! \brief Function that set the weight of the vertex
   *
   * \param id vertex id
   * \param weight to give to the vertex
   *
   */
  void setComputationCost(size_t id, size_t weight) {
    if (!verticesGotWeights) {
      verticesGotWeights = true;
    }

#ifdef SE_CLASS1  // question needed?
    if (id >= gp.getNVertex()) {
      std::cerr << __FILE__ << ":" << __LINE__
                << "Such vertex doesn't exist (id = " << id << ", "
                << "total size = " << gp.getNVertex() << ")\n";
    }
#endif

    // Update vertex in main graph
    gp.vertex(id).template get<nm_v_computation>() = weight;
  }

  /*! \brief Returns total number of neighbors of the sub-sub-domain id
   *
   * \param id id of the sub-sub-domain
   *
   * \return the number of neighborhood sub-sub-domains for each sub-domain
   *
   */
  size_t getNSubSubDomainNeighbors(size_t id) {
#ifdef SE_CLASS1
    if (id >= gp.getNVertex())
      std::cerr << __FILE__ << ":" << __LINE__
                << "Such vertex doesn't exist (id = " << id << ", "
                << "total size = " << gp.getNVertex() << ")\n";
#endif

    return gp.getNChilds(id);
  }

  /*! \brief Set the tolerance for each partition
   *
   * \param tol tolerance
   *
   */
  template <typename Graph>
  void setDistTol(Graph& graph, double tol) {
    graph.setDistTol(tol);
  }

  void setMigrationCosts(const float migration,
                         const size_t norm,
                         const size_t ts) {
    for (auto i = 0; i < getNSubSubDomains(); i++) {
      setMigrationCost(i, norm * migration);

      for (auto s = 0; s < getNSubSubDomainNeighbors(i); s++) {
        // We have to remove getSubSubDomainComputationCost(i) otherwise the
        // graph is not directed
        setCommunicationCost(i, s, 1 * ts);
      }
    }
  }

  template <typename SubDomains>
  void onEnd(SubDomains& sub_domains) {
    domain_nn_calculator_cart<dim>::reset();
    domain_nn_calculator_cart<dim>::setParameters(proc_box);
  }

  /*! \brief Get the current graph (main)
   *
   */
  Graph_CSR<nm_v<dim>, nm_e>& getGraph() { return gp; }

  template <typename Graph>
  void reset(Graph& graph) {
    if (is_distributed) {
      graph.reset(gp, vtxdist, m2g, verticesGotWeights);
    } else {
      graph.initSubGraph(gp, vtxdist, m2g, verticesGotWeights);
    }
  }

  /*! \brief Refine current decomposition
   *
   * It makes a refinement of the current decomposition using Parmetis function
   * RefineKWay After that it also does the remapping of the graph
   *
   */
  template <typename Graph>
  void refine(Graph& graph) {
    graph.reset(gp, vtxdist, m2g, verticesGotWeights);  // reset
    graph.refine(vtxdist);                              // refine
    distribute(graph);
  }

  /*! \brief Print the current distribution and save it to VTK file
   *
   * \param file filename
   *
   */
  void write(const std::string& file) {
    // todo VTKWriter<Graph_CSR<nm_v<dim>, nm_e>, VTK_GRAPH> gv2(gp);
    // todo gv2.write(std::to_string(v_cl.getProcessUnitID()) + "_" + file +
    // ".vtk");
  }

  openfpm::vector<rid>& getVtxdist() { return vtxdist; }

  /*! \brief Callback of the sendrecv to set the size of the array received
   *
   * \param msg_i Index of the message
   * \param total_msg Total numeber of messages
   * \param total_p Total number of processors to comunicate with
   * \param i Processor id
   * \param ri Request id
   * \param ptr Void pointer parameter for additional data to pass to the
   * call-back
   */
  static void* message_receive(size_t msg_i,
                               size_t total_msg,
                               size_t total_p,
                               size_t i,
                               size_t ri,
                               size_t tag,
                               void* ptr) {
    openfpm::vector<openfpm::vector<idx_t>>* v =
        static_cast<openfpm::vector<openfpm::vector<idx_t>>*>(ptr);

    v->get(i).resize(msg_i / sizeof(idx_t));

    return &(v->get(i).get(0));
  }

  /*! \brief Get the global id of the vertex given the re-mapped one
   *
   * \param remapped id
   * \return global id
   *
   */
  gid getVertexGlobalId(rid n) { return m2g.find(n)->second; }

  /*! \brief It update the full decomposition
   *
   *
   */
  template <typename Graph>
  void distribute(Graph& graph) {
    reset(graph);

    //! Get the processor id
    size_t p_id = v_cl.getProcessUnitID();

    //! Get the number of processing units
    size_t Np = v_cl.getProcessingUnits();

    // Number of local vertex
    size_t nl_vertex = vtxdist.get(p_id + 1).id - vtxdist.get(p_id).id;

    //! Get result partition for this processors
    idx_t* partition = graph.getPartition();

    //! Prepare vector of arrays to contain all partitions
    partitions.get(p_id).resize(nl_vertex);

    if (nl_vertex != 0) {
      std::copy(partition, partition + nl_vertex, &partitions.get(p_id).get(0));
    }

    // Reset data structure to keep trace of new vertices distribution in
    // processors (needed to update main graph)
    for (size_t i = 0; i < Np; ++i) {
      v_per_proc.get(i).clear();
    }

    // Communicate the local distribution to the other processors
    // to reconstruct individually the global graph
    openfpm::vector<size_t> prc;
    openfpm::vector<size_t> sz;
    openfpm::vector<void*> ptr;

    for (size_t i = 0; i < Np; i++) {
      if (i != v_cl.getProcessUnitID()) {
        partitions.get(i).clear();
        prc.add(i);
        sz.add(nl_vertex * sizeof(idx_t));
        ptr.add(partitions.get(p_id).getPointer());
      }
    }

    if (prc.size() == 0) {
      v_cl.sendrecvMultipleMessagesNBX(
          0, NULL, NULL, NULL, message_receive, &partitions, NONE);
    } else {
      v_cl.sendrecvMultipleMessagesNBX(prc.size(),
                                       &sz.get(0),
                                       &prc.get(0),
                                       &ptr.get(0),
                                       message_receive,
                                       &partitions,
                                       NONE);
    }

    // Update graphs with the received data
    updateGraphs();

    is_distributed = true;
  }

  /*! \brief Return the ghost
   *
   *
   * \return the ghost extension
   *
   */
  Ghost<dim, T>& getGhost() { return ghost; }

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

  /*! \brief Distribution grid
   *
   * \return the grid
   */
  DGrid getGrid() { return gr; }

  void setParameters(
      DGrid& grid_dec,
      const grid_sm<dim, void>& sec_dist = grid_sm<dim, void>()) {
    if (sec_dist.size(0) != 0) {
      gr.setDimensions(sec_dist.getSize());
    } else {
      gr = grid_dec;
    }
  }

  /*! \brief Create the Cartesian graph
   *
   * \param grid info
   * \param dom domain
   */
  void createCartGraph(const size_t (&bc)[dim], ::Box<dim, T>& domain) {
    // Create a cartesian grid graph
    CartesianGraphFactory<dim, Graph_CSR<nm_v<dim>, nm_e>> g_factory_part;
    gp = g_factory_part.template construct<NO_EDGE, nm_v_id, T, dim - 1, 0>(
        gr.getSize(), domain, bc);
    initLocalToGlobalMap();

    //! Get the number of processing units
    size_t Np = v_cl.getProcessingUnits();

    //! Division of vertices in Np graphs
    //! Put (div+1) vertices in mod graphs
    //! Put div vertices in the rest of the graphs
    size_t mod_v = gr.size() % Np;
    size_t div_v = gr.size() / Np;

    for (size_t i = 0; i <= Np; i++) {
      if (i < mod_v) {
        vtxdist.get(i).id = (div_v + 1) * i;
      } else {
        vtxdist.get(i).id = (div_v)*i + mod_v;
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

private:
  bool is_distributed = false;

  //! Processor domain bounding box
  ::Box<dim, size_t> proc_box;

  //! Global sub-sub-domain graph
  Graph_CSR<nm_v<dim>, nm_e> gp;

  //! ghost info
  Ghost<dim, T> ghost;

  //! Structure that store the cartesian grid information
  DGrid gr;

  //! Init vtxdist needed for Parmetis
  //
  // vtxdist is a common array across processor, it indicate how
  // vertex are distributed across processors
  //
  // Example we have 3 processors
  //
  // processor 0 has 3 vertices
  // processor 1 has 5 vertices
  // processor 2 has 4 vertices
  //
  // vtxdist contain, 0,3,8,12
  //
  // vtx dist is the unique global-id of the vertices
  //
  openfpm::vector<rid> vtxdist;

  //! partitions
  openfpm::vector<openfpm::vector<idx_t>> partitions;

  //! Init data structure to keep trace of new vertices distribution in
  //! processors (needed to update main graph)
  openfpm::vector<openfpm::vector<gid>> v_per_proc;

  //! Hashmap to access to the global position given the re-mapped one (needed
  //! for access the map)
  std::unordered_map<rid, gid> m2g;

  //! Id of the sub-sub-domain where we set the costs
  openfpm::vector<size_t> sub_sub_owner;

  //! Flag to check if weights are used on vertices
  bool verticesGotWeights = false;

  /*! \brief operator to remap vertex to a new position
   *
   * \param n re-mapped position
   * \param g global position
   *
   */
  void setMapId(rid n, gid g) { m2g[n] = g; }

  /*! \brief Update main graph ad subgraph with the received data of the
   * partitions from the other processors
   *
   */
  void updateGraphs() {
    sub_sub_owner.clear();

    size_t Np = v_cl.getProcessingUnits();

    // Init n_vtxdist to gather informations about the new decomposition
    openfpm::vector<rid> n_vtxdist(Np + 1);
    for (size_t i = 0; i <= Np; i++)
      n_vtxdist.get(i).id = 0;

    // Update the main graph with received data from processor i
    for (size_t i = 0; i < Np; i++) {
      size_t ndata = partitions.get(i).size();
      size_t k = 0;

      // Update the main graph with the received informations
      for (rid l = vtxdist.get(i); k < ndata && l < vtxdist.get(i + 1);
           k++, ++l) {
        // Create new n_vtxdist (just count processors vertices)
        ++n_vtxdist.get(partitions.get(i).get(k) + 1);

        // vertex id from vtx to grobal id
        auto v_id = m2g.find(l)->second.id;

        // Update proc id in the vertex (using the old map)
        gp.template vertex_p<nm_v_proc_id>(v_id) = partitions.get(i).get(k);

        if (partitions.get(i).get(k) == (long int)v_cl.getProcessUnitID()) {
          sub_sub_owner.add(v_id);
        }

        // Add vertex to temporary structure of distribution (needed to update
        // main graph)
        v_per_proc.get(partitions.get(i).get(k)).add(getVertexGlobalId(l));
      }
    }

    // Create new n_vtxdist (accumulate the counters)
    for (size_t i = 2; i <= Np; i++) {
      n_vtxdist.get(i) += n_vtxdist.get(i - 1);
    }

    // Copy the new decomposition in the main vtxdist
    for (size_t i = 0; i <= Np; i++) {
      vtxdist.get(i) = n_vtxdist.get(i);
    }

    openfpm::vector<size_t> cnt;
    cnt.resize(Np);

    for (size_t i = 0; i < gp.getNVertex(); ++i) {
      size_t pid = gp.template vertex_p<nm_v_proc_id>(i);

      rid j = rid(vtxdist.get(pid).id + cnt.get(pid));
      gid gi = gid(i);

      gp.template vertex_p<nm_v_id>(i) = j.id;
      cnt.get(pid)++;

      setMapId(j, gi);
    }
  }
};

#endif  // SRC_DECOMPOSITION_ABSTRACTDISTRIBUTIONSTRATEGY_HPP
