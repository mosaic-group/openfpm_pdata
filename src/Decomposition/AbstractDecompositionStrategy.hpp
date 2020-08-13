#ifndef SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP

#include <cmath>
#include <initializer_list>
#include <unordered_map>
#include <utility>
#include <vector>
#include "Decomposition/dec_optimizer.hpp"

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
#include "util/mathutil.hpp"
#include "util/se_util.hpp"
#include "util/generic.hpp"

template <unsigned int dim,
          typename T,
          typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin,
          typename DGrid = grid_sm<dim, void>>
class AbstractDecompositionStrategy
  : public ie_loc_ghost<dim, T, layout_base, Memory>,
    public nn_prcs<dim, T, layout_base, Memory>,
    public ie_ghost<dim, T, Memory, layout_base>,
    public domain_icell_calculator<dim, T, layout_base, Memory> {
  //! Type of the domain we are going to decompose
  using domain_type = T;

  //! It simplify to access the SpaceBox element
  using Box = SpaceBox<dim, T>;

  using SubDomains = openfpm::
      vector<Box, Memory, typename layout_base<Box>::type, layout_base>;

public:
  //! Structure that decompose the space into cells without creating them
  //! useful to convert positions to CellId or sub-domain id in this case
  CellDecomposer_sm<dim, T, shift<dim, T>> cd;  // todo private

  /*! \brief Abstract decomposition constructor
   *
   * \param v_cl Virtual cluster, used internally to handle or pipeline
   * comm0unication
   *
   */
  AbstractDecompositionStrategy(Vcluster<>& v_cl)
    : nn_prcs<dim, T, layout_base, Memory>(v_cl), v_cl(v_cl), ref_cnt(0) {
    bbox.zero();  // Reset the box to zero
  }

  //! Destructor
  ~AbstractDecompositionStrategy() {
    // question ref counter?
  }

  //! Increment the reference counter
  void incRef() { ref_cnt++; }

  //! Decrement the reference counter
  void decRef() { ref_cnt--; }

  //! Return the reference counter
  long int ref() { return ref_cnt; }

  /*! \brief Decomposition grid
   *
   * \return the grid
   */
  DGrid& getGrid() { return gr; }

  /*! \brief Calculate communication and migration costs
   */
  template <typename Ghost>
  std::pair<float, size_t> computeCommunicationCosts(Ghost& ghost) {
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

  /*! \brief Return the box of the physical domain
   *
   * \return The physical domain box
   *
   */
  const ::Box<dim, T>& getDomain() const { return domain; }

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
  ::Box<dim, T>& getProcessorBounds() { return bbox; }

  /*! \brief Start decomposition
   *
   */
  template <typename Model, typename Graph>
  void decompose(Model m, Graph& graph, openfpm::vector<rid>& vtxdist) {
    graph.decompose(vtxdist);  // decompose
  }

  template <typename Graph, typename Ghost>
  void merge(Graph& graph, Ghost& ghost, DGrid gr_dist) {
    createSubdomains(graph, ghost, gr_dist);
    calculateGhostBoxes(ghost);
  }

  template <typename Ghost>
  void onEnd(Ghost& ghost) {
    domain_icell_calculator<dim, T, layout_base, Memory>::
        CalculateInternalCells(
            v_cl,
            ie_ghost<dim, T, Memory, layout_base>::private_get_vb_int_box(),
            sub_domains,
            this->getProcessorBounds(),
            ghost.getRcut(),
            ghost);
  }

  bool shouldSetCosts() { return !costBeenSet; }

  SubDomains& getSubDomains() {
    return sub_domains;
  }

  void setParameters(
      const size_t (&div_)[dim],
      ::Box<dim, T>& domain_,
      const size_t (&bc)[dim],
      const grid_sm<dim, void>& sec_dist = grid_sm<dim, void>()) {
    for (size_t i = 0; i < dim; i++) {
      this->bc[i] = bc[i];  // todo std::copy
    }

    // Set the decomposition parameters
    gr.setDimensions(div_);
    domain = domain_;
    cd.setDimensions(domain, div_, 0);

    // calc magnification factor dec-dist
    calculate_magn(sec_dist);
  }
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

private:
  //! Structure that store the cartesian grid information
  DGrid gr;

  //! Box Spacing
  T spacing[dim];

  //! Magnification factor between distribution and
  //! decomposition
  size_t magn[dim];

  //! rectangular domain to decompose
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
  Vcluster<>& v_cl;

  //! the set of all local sub-domain as vector
  SubDomains sub_domains;

  //! the remote set of all sub-domains as vector of 'sub_domains' vectors
  mutable openfpm::vector<Box_map<dim, T>,
                          Memory,
                          typename layout_base<Box_map<dim, T>>::type,
                          layout_base>
      sub_domains_global;

  //! for each sub-domain, contain the list of the neighborhood processors
  openfpm::vector<openfpm::vector<long unsigned int>> box_nn_processor;

  //! Structure that contain for each sub-sub-domain box the processor id
  //! exist for efficient global communication
  CellList<dim, T, Mem_fast<Memory, int>, shift<dim, T>>
      fine_s;  // todo cellist to find particle

  //! set of Boxes produced by the decomposition optimizer
  openfpm::vector<::Box<dim, size_t>> loc_box;

  //! Processor domain bounding box
  ::Box<dim, size_t> proc_box;

  void initialize_fine_s(const ::Box<dim, T>& domain) {
    fine_s.clear();
    size_t div_g[dim];

    // We reduce the size of the cells by a factor 8 in 3d 4 in 2d
    for (size_t i = 0; i < dim; i++) {
      div_g[i] = (gr.size(i) == 1) ? 1 : gr.size(i) / 2;
    }

    fine_s.Initialize(domain, div_g);
  }

  /*! \brief It convert the box from the domain decomposition into sub-domain
   *
   * The decomposition box from the domain-decomposition contain the box in
   * integer coordinates. This box is converted into a continuos box. It also
   * adjust loc_box if the distribution grid and the decomposition grid are
   * different.
   *
   * \param loc_box local box
   *
   * \return the corresponding sub-domain
   *
   */
  template <typename Memory_bx>
  Box convertDecBoxIntoSubDomain(
      encapc<1, ::Box<dim, size_t>, Memory_bx> loc_box) {
    // A point with all coordinate to one
    size_t one[dim];
    for (size_t i = 0; i < dim; i++) {
      one[i] = 1;
    }

    SpaceBox<dim, size_t> sub_dc = loc_box;
    SpaceBox<dim, size_t> sub_dce = sub_dc;
    sub_dce.expand(one);
    sub_dce.mul(magn);

    // shrink by one
    for (size_t i = 0; i < dim; i++) {
      loc_box.template get<Box::p1>()[i] = sub_dce.getLow(i);
      loc_box.template get<Box::p2>()[i] = sub_dce.getHigh(i) - 1;
    }

    SpaceBox<dim, size_t> sub_d(sub_dce);
    sub_d.mul(spacing);
    sub_d += domain.getP1();

    // we add the

    // Fixing sub-domains to cover all the domain

    // Fixing sub_d
    // if (loc_box) is at the boundary we have to ensure that the box span the
    // full domain (avoiding rounding off error)
    for (size_t i = 0; i < dim; i++) {
      if (sub_dc.getHigh(i) == gr.size(i) - 1) {
        sub_d.setHigh(i, domain.getHigh(i));
      }

      if (sub_dc.getLow(i) == 0) {
        sub_d.setLow(i, domain.getLow(i));
      }
    }

    return sub_d;
  }

  void collect_all_sub_domains(
      openfpm::vector<Box_map<dim, T>,
                      Memory,
                      typename layout_base<Box_map<dim, T>>::type,
                      layout_base>& sub_domains_global) {
    sub_domains_global.clear();
    openfpm::vector<Box_map<dim, T>,
                    Memory,
                    typename layout_base<Box_map<dim, T>>::type,
                    layout_base>
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

  void construct_fine_s() {
    collect_all_sub_domains(sub_domains_global);

    // now draw all sub-domains in fine-s

    for (size_t i = 0; i < sub_domains_global.size(); i++) {
      // get the cells this box span
      const grid_key_dx<dim> p1 =
          fine_s.getCellGrid_me(sub_domains_global.template get<0>(i).getP1());
      const grid_key_dx<dim> p2 =
          fine_s.getCellGrid_pe(sub_domains_global.template get<0>(i).getP2());

      // Get the grid and the sub-iterator
      auto& gi = fine_s.getGrid();
      grid_key_dx_iterator_sub<dim> g_sub(gi, p1, p2);

      // add the box-id to the cell list
      while (g_sub.isNext()) {
        auto key = g_sub.get();
        fine_s.addCell(gi.LinId(key), i);

        ++g_sub;
      }
    }

    host_dev_transfer = false;
  }

  /*! \brief Calculate magnification
   *
   * \param gm distribution grid
   *
   */
  void calculate_magn(const grid_sm<dim, void>& gm) {
    if (gm.size() == 0) {
      for (size_t i = 0; i < dim; i++) {
        magn[i] = 1;
      }
    } else {
      for (size_t i = 0; i < dim; i++) {
        if (gr.size(i) % gm.size(i) != 0) {
          std::cerr << __FILE__ << ":" << __LINE__
                    << ".Error the decomposition grid specified as gr.size("
                    << i << ")=" << gr.size(i)
                    << " is not multiple of the distribution grid gm.size(" << i
                    << ")=" << gm.size(i) << std::endl;
        }

        magn[i] = gr.size(i) / gm.size(i);
      }
    }
  }

  /*! \brief Constructor, it decompose and distribute the sub-domains across the
   * processorscreateSubdomains
   *
   * \param v_cl Virtual cluster, used internally for communications
   * \param bc boundary conditions
   * \param opt option (one option is to construct)
   *
   */
  template <typename Graph, typename Ghost>
  void createSubdomains(Graph& graph,
                        Ghost& ghost,
                        DGrid gr_dist,
                        size_t opt = 0) {
    // Calculate the total number of box and and the spacing
    // on each direction
    // Get the box containing the domain
    SpaceBox<dim, T> bs = domain.getBox();

    for (unsigned int i = 0; i < dim; i++) {
      // Calculate the spacing
      spacing[i] = (bs.getHigh(i) - bs.getLow(i)) / gr.size(i);
    }

    // fill the structure that store the processor id for each sub-domain
    initialize_fine_s(domain);

    // Optimize the decomposition creating bigger spaces
    // And reducing Ghost over-stress
    dec_optimizer<dim, Graph_CSR<nm_v<dim>, nm_e>> d_o(graph,
                                                       gr_dist.getSize());

    // Ghost
    Ghost ghe;

    // Set the ghost
    printVar(dim);
    for (size_t i = 0; i < dim; i++) {
      ghe.setLow(i, static_cast<long int>(ghost.getLow(i) / spacing[i]) - 1);
      ghe.setHigh(i, static_cast<long int>(ghost.getHigh(i) / spacing[i]) + 1);
    }

    // optimize the decomposition or merge sub-sub-domain
    // breakpoint segfault
    d_o.template optimize<nm_v_sub_id, nm_v_proc_id>(
        graph, v_cl.getProcessUnitID(), loc_box, box_nn_processor, ghe, bc);

    // Initialize
    if (loc_box.size() > 0) {
      bbox = convertDecBoxIntoSubDomain(loc_box.get(0));
      proc_box = loc_box.get(0);
      sub_domains.add(bbox);
    } else {
      // invalidate all the boxes
      for (size_t i = 0; i < dim; i++) {
        proc_box.setLow(i, 0.0);
        proc_box.setHigh(i, 0);

        bbox.setLow(i, 0.0);
        bbox.setHigh(i, 0);
      }
    }

    // convert into sub-domain
    for (size_t s = 1; s < loc_box.size(); s++) {
      Box sub_d = convertDecBoxIntoSubDomain(loc_box.get(s));

      // add the sub-domain
      sub_domains.add(sub_d);

      // Calculate the bound box
      bbox.enclose(sub_d);
      proc_box.enclose(loc_box.get(s));
    }

    nn_prcs<dim, T, layout_base, Memory>::create(box_nn_processor, sub_domains);
    nn_prcs<dim, T, layout_base, Memory>::applyBC(domain, ghost, bc);

    // fill fine_s structure
    // fine_s structure contain the processor id for each sub-sub-domain
    // with sub-sub-domain we mean the sub-domain decomposition before
    // running dec_optimizer (before merging sub-domains)

    construct_fine_s();
    Initialize_geo_cell_lists();
  }

  void Initialize_geo_cell_lists() {
    // Get the processor bounding Box
    ::Box<dim, T> bound = getProcessorBounds();

    // Check if the box is valid
    if (bound.isValidN() == true) {
      // calculate the sub-divisions
      size_t div[dim];
      for (size_t i = 0; i < dim; i++) {
        div[i] = (size_t)((bound.getHigh(i) - bound.getLow(i)) /
                          cd.getCellBox().getP2()[i]);
      }

      // Initialize the geo_cell structure
      ie_ghost<dim, T, Memory, layout_base>::Initialize_geo_cell(bound, div);

      // Initialize shift vectors
      ie_ghost<dim, T, Memory, layout_base>::generateShiftVectors(domain, bc);
    }
  }

  /*! \brief It calculate the internal ghost boxes
   */

  template <typename Ghost>
  void calculateGhostBoxes(Ghost& ghost) {
    // Intersect all the local sub-domains with the sub-domains of the
    // contiguous processors

    // create the internal structures that store ghost information
    ie_ghost<dim, T, Memory, layout_base>::create_box_nn_processor_ext(
        v_cl, ghost, sub_domains, box_nn_processor, *this);
    ie_ghost<dim, T, Memory, layout_base>::create_box_nn_processor_int(
        v_cl, ghost, sub_domains, box_nn_processor, *this);

    ie_loc_ghost<dim, T, layout_base, Memory>::create(
        sub_domains, domain, ghost, bc);

    // Ghost box information must be re-offloaded
    host_dev_transfer = false;
    ie_ghost<dim, T, Memory, layout_base>::reset_host_dev_transfer();
  }
};
#endif  // SRC_DECOMPOSITION_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
