#ifndef SRC_DECOMPOSITION_ABSTRACTDISTRIBUTIONSTRATEGY_HPP
#define SRC_DECOMPOSITION_ABSTRACTDISTRIBUTIONSTRATEGY_HPP

#include "Graph/CartesianGraphFactory.hpp"
#include "Graph/ids.hpp"
#include "SubdomainGraphNodes.hpp"

/*! \brief Class that distribute sub-sub-domains across processors
 */
template <unsigned int dim, typename T, typename DistributionGraph, typename Memory = HeapMemory, template <typename> class layout_base = memory_traits_lin>
class AbstractDistributionStrategy {

public:
  //! Vcluster
  Vcluster<> &v_cl;

  /*! Constructor
   *
   * \param v_cl Vcluster to use as communication object in this class
   */
  AbstractDistributionStrategy(Vcluster<> &v_cl)
      : is_distributed(false), v_cl(v_cl), graph(v_cl, v_cl.getProcessingUnits()),
        vtxdist(v_cl.getProcessingUnits() + 1),
        _partitions(v_cl.getProcessingUnits()),
        _v_per_proc(v_cl.getProcessingUnits()) {}

  /*! \brief Set the tolerance for each partition
   *
   * \param tol tolerance
   *
   */
  void setDistTol(double tol) {
    graph.setDistTol(tol);
  }

  /*! \brief It update the full decomposition
   */
  void distribute() {
    graph.decompose(getVtxdist());

    is_distributed = true;
  }

  void onEnd() {}

  template <typename DecompositionGraph>
  void reset(DecompositionGraph& gp) {
    if (is_distributed) {
      getGraph().reset(gp, getVtxdist(), m2g, verticesGotWeights);
    } else {
      getGraph().initSubGraph(gp, getVtxdist(), m2g, verticesGotWeights);
    }
  }

  template <typename DecompositionGraph>
  void refine(DecompositionGraph& gp) {
    graph.reset(gp, vtxdist, m2g, verticesGotWeights); // reset
    graph.refine(getVtxdist());
  }

  /*! \brief Print the current distribution and save it to VTK file
   *
   * \param file filename
   *
   */
  void write(const std::string &file) {
    // todo
  }

  /*! \brief operator to init ids vector
   *
   * operator to init ids vector
   *
   */
  // todo maybe DecompositionGraph as a template
  template <typename DecompositionGraph>
  void initLocalToGlobalMap(DecompositionGraph& gp) {
    gid g;
    rid i;
    i.id = 0;

    m2g.clear();
    for (; (size_t)i.id < gp.getNVertex(); ++i) {
      g.id = i.id;

      m2g.insert({i, g});
    }
  }

  openfpm::vector<rid> &getVtxdist() { return vtxdist; }

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
  static void *message_receive(size_t msg_i, size_t total_msg, size_t total_p,
                               size_t i, size_t ri, size_t tag, void *ptr) {
    openfpm::vector<openfpm::vector<idx_t>> *v =
        static_cast<openfpm::vector<openfpm::vector<idx_t>> *>(ptr);

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

  /*! \brief operator to remap vertex to a new position
   *
   * \param n re-mapped position
   * \param g global position
   *
   */
  void setMapId(rid n, gid g) { m2g[n] = g; }

  Vcluster<> & getVcluster() { return v_cl; }

  DistributionGraph & getGraph() { return graph; }

  openfpm::vector<openfpm::vector<idx_t>> & partitions() { return _partitions; }

  openfpm::vector<openfpm::vector<gid>> & v_per_proc() { return _v_per_proc; }

  // todo private
  //! Hashmap to access to the global position given the re-mapped one (needed
  //! for access the map)
  std::unordered_map<rid, gid> m2g;

private:
  bool is_distributed = false;

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
  openfpm::vector<openfpm::vector<idx_t>> _partitions;

  //! Init data structure to keep trace of new vertices distribution in
  //! processors (needed to update main graph)
  openfpm::vector<openfpm::vector<gid>> _v_per_proc;

  // graph of PUs (= processing unit(s))
  DistributionGraph graph;

  //! Flag to check if weights are used on vertices
  bool verticesGotWeights = false;
};

#endif // SRC_DECOMPOSITION_ABSTRACTDISTRIBUTIONSTRATEGY_HPP
