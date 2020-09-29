#ifndef SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP

#include <cmath>
#include <initializer_list>
#include <unordered_map>
#include <utility>
#include <vector>

#include "DLB/DLB.hpp"
#include "Decomposition/Decomposition.hpp"
#include "Decomposition/Domain_icells_cart.hpp"
#include "Decomposition/common.hpp"
#include "Decomposition/ie_ghost.hpp"
#include "Decomposition/ie_loc_ghost.hpp"
#include "Decomposition/nn_processor.hpp"
#include "Graph/CartesianGraphFactory.hpp"
#include "GraphMLWriter/GraphMLWriter.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include "NN/CellList/CellList.hpp"
#include "Space/Ghost.hpp"
#include "Space/Shape/Box.hpp"
#include "Space/Shape/Point.hpp"
#include "SubdomainGraphNodes.hpp"
#include "VCluster/VCluster.hpp"
#include "Vector/map_vector.hpp"
#include "config.h"
#include "data_type/aggregate.hpp"
#include "util/generic.hpp"
#include "util/mathutil.hpp"
#include "util/se_util.hpp"

template <unsigned int dim, typename T, typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin,
          typename DGrid = grid_sm<dim, void>>
class AbstractDecompositionStrategy
    : public ie_loc_ghost<dim, T, layout_base, Memory>,
      public nn_prcs<dim, T, layout_base, Memory>,
      public ie_ghost<dim, T, Memory, layout_base>,
      public domain_icell_calculator<dim, T, layout_base, Memory> {
  //! Type of the domain we are going to decompose
  using domain_type = T;
  using Box = SpaceBox<dim, T>;
  using SubDomains =
      openfpm::vector<Box, Memory, typename layout_base<Box>::type,
                      layout_base>;

public:
  /*! \brief Abstract decomposition constructor
   *
   * \param v_cl Virtual cluster, used internally to handle or pipeline
   * comm0unication
   *
   */
  AbstractDecompositionStrategy(Vcluster<> &v_cl)
      : nn_prcs<dim, T, layout_base, Memory>(v_cl), v_cl(v_cl), ref_cnt(0) {
    bbox.zero(); // Reset the box to zero
  }

  /*! \brief Copy constructor
   *
   * \param cart object to copy
   *
   */
  AbstractDecompositionStrategy(
      const AbstractDecompositionStrategy<dim, T, Memory, layout_base, DGrid>
          &cart)
      : nn_prcs<dim, T, layout_base, Memory>(cart.v_cl), v_cl(cart.v_cl),
        ref_cnt(0) {
    this->operator=(cart);
  }

  /*! \brief Copy the element
   *
   * \param cart element to copy
   *
   * \return itself
   *
   */
  AbstractDecompositionStrategy<dim, T, Memory, layout_base, DGrid> &
  operator=(const AbstractDecompositionStrategy &cart) {
    static_cast<ie_loc_ghost<dim, T, layout_base, Memory> *>(this)->operator=(
        static_cast<ie_loc_ghost<dim, T, layout_base, Memory>>(cart));
    static_cast<nn_prcs<dim, T, layout_base, Memory> *>(this)->operator=(
        static_cast<nn_prcs<dim, T, layout_base, Memory>>(cart));
    static_cast<ie_ghost<dim, T, Memory, layout_base> *>(this)->operator=(
        static_cast<ie_ghost<dim, T, Memory, layout_base>>(cart));
    sub_domains = cart.sub_domains;
    box_nn_processor = cart.box_nn_processor;
    fine_s = cart.fine_s;
    cd = cart.cd;
    domain = cart.domain;
    sub_domains_global = cart.sub_domains_global;

    for (size_t i = 0; i < dim; i++) {
      spacing[i] = cart.spacing[i];
      magn[i] = cart.magn[i];
    };

    bbox = cart.bbox;

    for (size_t i = 0; i < dim; i++) {
      bc[i] = cart.bc[i];
    }

    return *this;
  }

  //! Destructor
  ~AbstractDecompositionStrategy() {}

  /*! \brief class to select the returned id by ghost_processorID
   *
   */
  class box_id {
  public:
    /*! \brief Return the box id
     *
     * \param p structure containing the id informations
     * \param b_id box_id
     *
     * \return box id
     *
     */
    inline static size_t id(p_box<dim, T> &p, size_t b_id) { return b_id; }
  };

  /*! \brief class to select the returned id by ghost_processorID
   *
   */
  class processor_id {
  public:
    /*! \brief Return the processor id
     *
     * \param p structure containing the id informations
     * \param b_id box_id
     *
     * \return processor id
     *
     */
    template <typename encap_type>
    inline static size_t id(const encap_type &p, size_t b_id) {
      return p.template get<proc_>();
    }
  };

  /*! \brief class to select the returned id by ghost_processorID
   *
   */
  class lc_processor_id {
  public:
    /*! \brief Return the near processor id
     *
     * \param p structure containing the id informations
     * \param b_id box_id
     *
     * \return local processor id
     *
     */
    template <typename encap_type>
    inline static size_t id(const encap_type &p, size_t b_id) {
      return p.template get<lc_proc_>();
    }
  };

  /*! \brief class to select the returned id by ghost_processorID
   *
   */
  class shift_id {
  public:
    /*! \brief Return the shift id
     *
     * \param p structure containing the id informations
     * \param b_id box_id
     *
     * \return shift_id id
     *
     */
    template <typename encap_type>
    inline static size_t id(const encap_type &p, size_t b_id) {
      return p.template get<shift_id_>();
    }
  };

  //! Increment the reference counter
  void incRef() { ref_cnt++; }

  //! Decrement the reference counter
  void decRef() { ref_cnt--; }

  //! Return the reference counter
  long int ref() { return ref_cnt; }

  /*! \brief Calculate communication and migration costs
   */
  template <typename Ghost>
  std::pair<float, size_t> computeCommunicationCosts(Ghost &ghost) {
    const Box cellBox = cd.getCellBox();
    const float b_s = static_cast<float>(cellBox.getHigh(0));
    const float gh_s = static_cast<float>(ghost.getHigh(0));

    // compute the gh_area for 2 dim case
    float gh_v = (gh_s * b_s);

    // multiply for sub-sub-domain side for each domain
    for (auto i = 2; i < dim; i++) {
      gh_v *= b_s;
    }

    const size_t norm = (size_t)(1.0 / gh_v);
    float migration = pow(b_s, dim);

    costBeenSet = true;

    return std::make_pair(migration, norm);
  }

  bool shouldSetCosts() { return !costBeenSet; }

  /*! \brief Return the box of the physical domain
   *
   * \return The physical domain box
   *
   */
  const ::Box<dim, T> &getDomain() const { return domain; }

  /*! \brief Delete the decomposition and reset the data-structure
   *
   *
   */
  void reset() {
    sub_domains.clear();
    box_nn_processor.clear();
    fine_s.clear();
    loc_box.clear();
    nn_prcs<dim, T, layout_base, Memory>::reset();
    ie_ghost<dim, T, Memory, layout_base>::reset();
    ie_loc_ghost<dim, T, layout_base, Memory>::reset();
  }

  /*! \brief Return the bounding box containing union of all the sub-domains for
   * the local processor
   *
   * \return The bounding box
   *
   */
  ::Box<dim, T> &getProcessorBounds() { return bbox; }

  template <typename Model, typename Graph>
  void decompose(Model m, Graph &graph, openfpm::vector<rid> &vtxdist) {
    graph.decompose(vtxdist); // decompose
  }

  template <typename Graph>
  void merge(Graph &graph, Ghost<dim, T> &ghost, DGrid gr_dist) {}

  void onEnd(Ghost<dim, T> &ghost) {
    domain_icell_calculator<dim, T, layout_base, Memory>::
        CalculateInternalCells(
            v_cl,
            ie_ghost<dim, T, Memory, layout_base>::private_get_vb_int_box(),
            sub_domains, this->getProcessorBounds(), ghost.getRcut(), ghost);
  }

  SubDomains &getSubDomains() { return sub_domains; }

  /*! \brief Get the number of local sub-domains
   *
   * \return the number of sub-domains
   *
   */
  size_t getNSubDomain() { return sub_domains.size(); }

  /*! \brief function to check the consistency of the information of the
   * decomposition
   *
   * \return false if is inconsistent
   *
   */
  bool check_consistency() {
    return ie_loc_ghost<dim, T, layout_base, Memory>::check_consistency(
        getNSubDomain());
  }

// todo private:
  //! Structure that decompose the space into cells without creating them
  //! useful to convert positions to CellId or sub-domain id in this case
  CellDecomposer_sm<dim, T, shift<dim, T>> cd;  // todo abstract ???

  //! Box Spacing
  T spacing[dim];

  //! Magnification factor between distribution and
  //! decomposition
  size_t magn[dim];

  //! (rectangular) domain to decompose
  ::Box<dim, T> domain;

  //! Processor bounding box
  ::Box<dim, T> bbox;

  //! reference counter of the object in case is shared between object
  long int ref_cnt;

  //! bool that indicate whenever the buffer has been already transfer to device
  bool host_dev_transfer = false;

  //! Boundary condition info
  size_t bc[dim];

  bool costBeenSet = false;

  //! Runtime virtual cluster machine
  Vcluster<> &v_cl;

  //! the set of all local sub-domain as vector
  SubDomains sub_domains;

  //! the remote set of all sub-domains as vector of 'sub_domains' vectors
  mutable openfpm::vector<Box_map<dim, T>, Memory,
                          typename layout_base<Box_map<dim, T>>::type,
                          layout_base>
      sub_domains_global;

  //! for each sub-domain, contain the list of the neighborhood processors
  openfpm::vector<openfpm::vector<long unsigned int>> box_nn_processor;

  //! Structure that contain for each sub-sub-domain box the processor id
  //! exist for efficient global communication
  CellList<dim, T, Mem_fast<Memory, int>, shift<dim, T>> fine_s;

  //! set of Boxes produced by the decomposition optimizer
  openfpm::vector<::Box<dim, size_t>> loc_box;

  //! Processor domain bounding box
  ::Box<dim, size_t> proc_box;

// todo private:
  void collect_all_sub_domains(
      openfpm::vector<Box_map<dim, T>, Memory,
                      typename layout_base<Box_map<dim, T>>::type, layout_base>
          &sub_domains_global) {
    sub_domains_global.clear();
    openfpm::vector<Box_map<dim, T>, Memory,
                    typename layout_base<Box_map<dim, T>>::type, layout_base>
        bm;

    for (size_t i = 0; i < sub_domains.size(); i++) {
      bm.add();

      bm.template get<0>(bm.size() - 1) =
          ::SpaceBox<dim, T>(sub_domains.get(i));
      bm.template get<1>(bm.size() - 1) = v_cl.rank();
    }

    v_cl.SGather<decltype(bm), decltype(sub_domains_global), layout_base>(
        bm, sub_domains_global, 0);

    size_t size = sub_domains_global.size();

    v_cl.max(size);
    v_cl.execute();

    sub_domains_global.resize(size);

    v_cl.Bcast(sub_domains_global, 0);
    v_cl.execute();
  }

  /*! \brief It calculate the internal ghost boxes
   */

  template <typename Ghost> void calculateGhostBoxes(Ghost &ghost) {
    // Intersect all the local sub-domains with the sub-domains of the
    // contiguous processors

    // create the internal structures that store ghost information
    ie_ghost<dim, T, Memory, layout_base>::create_box_nn_processor_ext(
        v_cl, ghost, sub_domains, box_nn_processor, *this);
    ie_ghost<dim, T, Memory, layout_base>::create_box_nn_processor_int(
        v_cl, ghost, sub_domains, box_nn_processor, *this);

    ie_loc_ghost<dim, T, layout_base, Memory>::create(sub_domains, domain,
                                                      ghost, bc);

    // Ghost box information must be re-offloaded
    host_dev_transfer = false;
    ie_ghost<dim, T, Memory, layout_base>::reset_host_dev_transfer();
  }
};
#endif // SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP