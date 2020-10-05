#ifndef SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP

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
#include "util/generic.hpp"
#include "util/mathutil.hpp"
#include "util/se_util.hpp"

#include "Decomposition/AbstractDecompositionStrategy.hpp"

template <unsigned int dim, typename domain_type, typename Memory = HeapMemory, template <typename> class layout_base = memory_traits_lin,
          typename AbstractDecStrategy =
              AbstractDecompositionStrategy<dim, domain_type>>
class CartDecompositionStrategy
    : public ie_loc_ghost<dim, domain_type, layout_base, Memory>,
      public nn_prcs<dim, domain_type, layout_base, Memory>,
      public ie_ghost<dim, domain_type, Memory, layout_base>,
      public domain_icell_calculator<dim, domain_type, layout_base, Memory> {
  // todo reuse AbstractDecompositionStrategy types
  using Box = SpaceBox<dim, domain_type>;
  using SubDomains =
      openfpm::vector<Box, Memory, typename layout_base<Box>::type,
                      layout_base>;

public:
    //! Helpers
  /*! \brief class to select the returned id by ghost_processorID
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
    inline static size_t id(p_box<dim, domain_type> &p, size_t b_id) { return b_id; }
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

  CartDecompositionStrategy(Vcluster<> &v_cl)
      : nn_prcs<dim, domain_type, layout_base, Memory>(v_cl), _inner(v_cl) {
        bbox.zero(); // Reset the box to zero
  }

  /*! \brief Copy constructor
   *
   * \param cart object to copy
   *
   */
  CartDecompositionStrategy(
      const CartDecompositionStrategy<dim, domain_type, AbstractDecStrategy>& cart) : nn_prcs<dim, domain_type, layout_base, Memory>(cart.getVcluster()) {
    this->operator=(cart);
  }

  /*! \brief Copy the element
   *
   * \param cart element to copy
   *
   * \return itself
   *
   */
  CartDecompositionStrategy<dim, domain_type> &
  operator=(const CartDecompositionStrategy &cart) {
    // todo
  }

  //! Destructor
  ~CartDecompositionStrategy() {}

  /*! \brief Return the bounding box containing union of all the sub-domains for
   * the local processor
   *
   * \return The bounding box
   *
   */
  ::Box<dim, domain_type> &getProcessorBounds() { return bbox; }

  /*! \brief Calculate magnification
   *
   * \param gm distribution grid
   *
   */
  void calculate_magn(const grid_sm<dim,void> & gm) {
    if (gm.size() == 0)
    {
      for (size_t i = 0 ; i < dim ; i++)
        magn[i] = 1;
    }
    else
    {
      for (size_t i = 0 ; i < dim ; i++)
      {
        if (gr.size(i) % gm.size(i) != 0)
          std::cerr << __FILE__ << ":" << __LINE__ << ".Error the decomposition grid specified as gr.size(" << i << ")=" << gr.size(i) << " is not multiple of the distribution grid gm.size(" << i << ")=" << gm.size(i) << std::endl;

        magn[i] = gr.size(i) / gm.size(i);
      }
    }
  }

  void setParameters(
    const size_t (&div_)[dim],
    ::Box<dim, domain_type> &domain_,
    const size_t (&bc)[dim],
    const Ghost<dim, domain_type> &ghost,
    const grid_sm<dim, void> &sec_dist = grid_sm<dim, void>()) {
    inner().setParameters(domain_, bc, ghost);
    gr.setDimensions(div_);
    cd.setDimensions(inner().getDomain(), div_, 0);

    // calc magnification factor dec-dist
    calculate_magn(sec_dist);
  }

  std::pair<float, size_t> computeCommunicationCosts() {
    const Box cellBox = cd.getCellBox();
    const float b_s = static_cast<float>(cellBox.getHigh(0));
    const float gh_s = static_cast<float>(inner().getGhost().getHigh(0));

    // compute the gh_area for 2 dim case
    float gh_v = (gh_s * b_s);

    // multiply for sub-sub-domain side for each domain
    for (auto i = 2; i < dim; i++) {
      gh_v *= b_s;
    }

    const size_t norm = (size_t)(1.0 / gh_v);
    float migration = pow(b_s, dim);

    inner().computeCommunicationCosts();

    return std::make_pair(migration, norm);
  }

  void decompose() {
    reset();

    if (inner().shouldSetCosts()) {
      computeCommunicationAndMigrationCosts(1);
    }
  }

  void onEnd() {
    domain_icell_calculator<dim, domain_type, layout_base, Memory>::
        CalculateInternalCells(
            inner().getVcluster(),
            ie_ghost<dim, domain_type, Memory, layout_base>::private_get_vb_int_box(),
            sub_domains, this->getProcessorBounds(), inner().getGhost().getRcut(), inner().getGhost());
  }

  void reset() {
    sub_domains.clear();
    box_nn_processor.clear();
    fine_s.clear();
    loc_box.clear();

    nn_prcs<dim, domain_type, layout_base, Memory>::reset();
    ie_ghost<dim, domain_type, Memory, layout_base>::reset();
    ie_loc_ghost<dim, domain_type, layout_base, Memory>::reset();

    inner().reset();
  }

  void computeCommunicationAndMigrationCosts(size_t ts) {
    float migration = 0;

    SpaceBox<dim, domain_type> cellBox = cd.getCellBox();
    float b_s = static_cast<float>(cellBox.getHigh(0));
    float gh_s = static_cast<float>(inner().getGhost().getHigh(0));

    // compute the gh_area for 2 dim case
    float gh_v = (gh_s * b_s);

    // multiply for sub-sub-domain side for each domain
    for (size_t i = 2; i < dim; i++)
    {
      /* coverity[dead_error_line] */
      gh_v *= b_s;
    }

    size_t norm = (size_t) (1.0 / gh_v);

    migration = pow(b_s, dim);

    size_t prev = 0;

    for (size_t i = 0; i < getNSubSubDomains(); i++)
    {
      setMigrationCost(i, norm * migration /* * getSubSubDomainComputationCost(i)*/ );

      for (size_t s = 0; s < getNSubSubDomainNeighbors(i); s++)
      {
        // We have to remove getSubSubDomainComputationCost(i) otherwise the graph is
        // not directed
        setCommunicationCost(i, s, 1 /** getSubSubDomainComputationCost(i)*/  *  ts);
      }
      prev += getNSubSubDomainNeighbors(i);
    }

    inner().computeCommunicationCosts();
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
    return ie_loc_ghost<dim, domain_type, layout_base, Memory>::check_consistency(
        getNSubDomain());
  }

  void initialize_fine_s(const ::Box<dim, domain_type> &domain) {
    fine_s.clear();
    size_t div_g[dim];

    // We reduce the size of the cells by a factor 8 in 3d 4 in 2d
    for (size_t i = 0; i < dim; ++i) {
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
  SpaceBox<dim, domain_type> convertDecBoxIntoSubDomain(
      encapc<1, ::Box<dim, size_t>, Memory_bx> loc_box) {
    // A point with all coordinate to one
    size_t one[dim];
    for (size_t i = 0; i < dim; i++) {
      one[i] = 1;  // todo use std::fill
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

    SpaceBox<dim, domain_type> sub_d(sub_dce);
    sub_d.mul(inner().spacing);
    sub_d += inner().getDomain().getP1();

    // we add the fixing sub-domains to cover all the domain
    // if (loc_box) is at the boundary we have to ensure that the box span the
    // full domain (avoiding rounding off error)
    for (size_t i = 0; i < dim; i++) {
      if (sub_dc.getHigh(i) == gr.size(i) - 1) {
        sub_d.setHigh(i, inner().getDomain().getHigh(i));
      }

      if (sub_dc.getLow(i) == 0) {
        sub_d.setLow(i, inner().getDomain().getLow(i));
      }
    }

    return sub_d;
  }

  /*! \brief Constructor, it decompose and distribute the sub-domains across the
   * processorscreateSubdomains
   *
   * \param v_cl Virtual cluster, used internally for communications
   * \param inner().bc boundary conditions
   * \param opt option (one option is to construct)
   *
   */
  void createSubdomains(grid_sm<dim, void> & gr_dist, size_t opt = 0) {
    // Calculate the total number of box and and the spacing
    // on each direction
    // Get the box containing the domain
    SpaceBox<dim, domain_type> bs = inner().getDomain().getBox();

    for (unsigned int i = 0; i < dim; i++) {
      // Calculate the spacing
      inner().spacing[i] = (bs.getHigh(i) - bs.getLow(i)) / gr.size(i);
    }

    // fill the structure that store the processor id for each sub-domain
    initialize_fine_s(inner().getDomain());

    // Optimize the decomposition creating bigger spaces
    // And reducing Ghost over-stress
    dec_optimizer<dim, Graph_CSR<nm_v<dim>, nm_e>> d_o(gp,
                                                       gr_dist.getSize());

    // Ghost
    Ghost<dim, long int> ghe;

    // Set the ghost
    for (size_t i = 0; i < dim; i++)
    {
      ghe.setLow(i, static_cast<long int>(inner().getGhost().getLow(i) / inner().spacing[i]) - 1);  // -1
      ghe.setHigh(i, static_cast<long int>(inner().getGhost().getHigh(i) / inner().spacing[i]) + 1);  // +1
    }

    // optimize the decomposition or merge sub-sub-domain
    d_o.template optimize<nm_v_sub_id, nm_v_proc_id>(
        gp, inner().getVcluster().getProcessUnitID(), loc_box, box_nn_processor, ghe, inner().bc);

    // Initialize
    if (loc_box.size() > 0)
    {
      bbox = convertDecBoxIntoSubDomain(loc_box.get(0));
      proc_box = loc_box.get(0);
      sub_domains.add(bbox);
    }
    else
    {
      // invalidate all the boxes
      for (size_t i = 0; i < dim; i++) {
        proc_box.setLow(i, 0.0);
        proc_box.setHigh(i, 0);

        bbox.setLow(i, 0.0);
        bbox.setHigh(i, 0);
      }
    }

    // convert into sub-domain
		for (size_t s = 1; s < loc_box.size(); s++)
		{
			SpaceBox<dim, domain_type> sub_d = convertDecBoxIntoSubDomain(loc_box.get(s));

			// add the sub-domain
			sub_domains.add(sub_d);

			// Calculate the bound box
			bbox.enclose(sub_d);
			proc_box.enclose(loc_box.get(s));
		}

    nn_prcs<dim, domain_type, layout_base, Memory>::create(box_nn_processor, sub_domains);
    nn_prcs<dim, domain_type, layout_base, Memory>::applyBC(inner().getDomain(), inner().getGhost(), inner().bc);

    // fill fine_s structure
    // fine_s structure contain the processor id for each sub-sub-domain
    // with sub-sub-domain we mean the sub-domain decomposition before
    // running dec_optimizer (before merging sub-domains)

    construct_fine_s();
    Initialize_geo_cell_lists();
  }

  void construct_fine_s() {
    collect_all_sub_domains(sub_domains_global);

    // now draw all sub-domains in fine-s
    for (size_t i = 0; i < sub_domains_global.size(); ++i) {
      // get the cells this box span
      const grid_key_dx<dim> p1 =
          fine_s.getCellGrid_me(sub_domains_global.template get<0>(i).getP1());
      const grid_key_dx<dim> p2 =
          fine_s.getCellGrid_pe(sub_domains_global.template get<0>(i).getP2());

      // Get the grid and the sub-iterator
      auto &gi = fine_s.getGrid();
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

  void Initialize_geo_cell_lists() {
    // Get the processor bounding Box
    ::Box<dim, domain_type> bound = getProcessorBounds();

    // Check if the box is valid
    if (bound.isValidN()) {
      // calculate the sub-divisions
      size_t div[dim];
      for (size_t i = 0; i < dim; ++i) {
        div[i] = (size_t)((bound.getHigh(i) - bound.getLow(i)) /
                          cd.getCellBox().getP2()[i]);
      }

      // Initialize the geo_cell structure
      ie_ghost<dim, domain_type, Memory, layout_base>::Initialize_geo_cell(bound, div);

      // Initialize shift vectors
      ie_ghost<dim, domain_type, Memory, layout_base>::generateShiftVectors(inner().getDomain(), inner().bc);
    }
  }

  void merge(grid_sm<dim, void>& gr_dist) {
    createSubdomains(gr_dist);
    calculateGhostBoxes();
  }

  void collect_all_sub_domains(
      openfpm::vector<Box_map<dim, domain_type>, Memory,
                      typename layout_base<Box_map<dim, domain_type>>::type, layout_base>
          &sub_domains_global) {
    sub_domains_global.clear();
    openfpm::vector<Box_map<dim, domain_type>, Memory,
                    typename layout_base<Box_map<dim, domain_type>>::type, layout_base>
        bm;

    for (size_t i = 0; i < sub_domains.size(); i++) {
      bm.add();

      bm.template get<0>(bm.size() - 1) =
          ::SpaceBox<dim, domain_type>(sub_domains.get(i));
      bm.template get<1>(bm.size() - 1) = inner().getVcluster().rank();
    }

    // todo inner().getVcluster().SGather<decltype(bm), decltype(sub_domains_global), layout_base>(
    // todo     bm, sub_domains_global, 0);

    size_t size = sub_domains_global.size();

    inner().getVcluster().max(size);
    inner().getVcluster().execute();

    sub_domains_global.resize(size);

    inner().getVcluster().Bcast(sub_domains_global, 0);
    inner().getVcluster().execute();
  }

  /*! \brief Set migration cost of the vertex id
   *
   * \param id of the vertex to update
   * \param migration cost of the migration
   */
  void setMigrationCost(size_t id, size_t migration) {
    gp.vertex(id).template get<nm_v_migration>() = migration;
  }

  void setMigrationCosts(const float migration, const size_t norm,
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

  /*! \brief Set communication cost of the edge id
   *
   * \param v_id Id of the source vertex of the edge
   * \param e i child of the vertex
   * \param communication Communication value
   */
  void setCommunicationCost(size_t v_id, size_t e, size_t communication) {
    gp.getChildEdge(v_id, e).template get<nm_e::communication>() =
        communication;
  }

  /*! \brief function that return the position of the vertex in the space
   *
   * \param id vertex id
   * \param pos vector that will contain x, y, z
   *
   */
  void getSubSubDomainPosition(size_t id, domain_type (&pos)[dim]) {
    // Copy the geometrical informations inside the pos vector
    pos[0] = gp.vertex(id).template get<nm_v_x>()[0];
    pos[1] = gp.vertex(id).template get<nm_v_x>()[1];
    if (dim == 3) {
      pos[2] = gp.vertex(id).template get<nm_v_x>()[2];
    }
  }

  void addComputationCost(size_t gid, size_t i) {
    size_t c = getSubSubDomainComputationCost(gid);
    setComputationCost(gid, c + i);
  }

  /*! \brief Function that set the weight of the vertex
   *
   * \param id vertex id
   * \param weight to give to the vertex
   *
   */
  void setComputationCost(size_t id, size_t weight) {
    // todo if (!verticesGotWeights) {
    // todo   verticesGotWeights = true;
    // todo }

    // Update vertex in main graph
    gp.vertex(id).template get<nm_v_computation>() = weight;
  }

  /*! \brief function that get the weight of the vertex
   * (computation cost of the sub-sub-domain id)
   *
   * \param id vertex id
   *
   */
  size_t getSubSubDomainComputationCost(size_t id) {
    return gp.vertex(id).template get<nm_v_computation>();
  }

  /*! \brief Returns total number of sub-sub-domains in the distribution graph
   *
   * \return the total number of sub-sub-domains
   *
   */
  size_t getNSubSubDomains() const { return gp.getNVertex(); }

  /*! \brief Add computation cost i to the subsubdomain with global id gid
   *
   * \param gid global id of the subsubdomain to update
   * \param i Cost increment
   */

  /*! \brief Returns total number of neighbors of the sub-sub-domain id
   *
   * \param id id of the sub-sub-domain
   *
   * \return the number of neighborhood sub-sub-domains for each sub-domain
   *
   */
  size_t getNSubSubDomainNeighbors(size_t id) {
    return gp.getNChilds(id);
  }

  /*! \brief It calculate the internal ghost boxes
   */

  // todo `use ghost ...`
  void calculateGhostBoxes() {
    // Intersect all the local sub-domains with the sub-domains of the
    // contiguous processors

    // create the internal structures that store ghost information
    ie_ghost<dim, domain_type, Memory, layout_base>::create_box_nn_processor_ext(
        inner().getVcluster(), inner().getGhost(), sub_domains, box_nn_processor, *this);
    ie_ghost<dim, domain_type, Memory, layout_base>::create_box_nn_processor_int(
        inner().getVcluster(), inner().getGhost(), sub_domains, box_nn_processor, *this);

    ie_loc_ghost<dim, domain_type, layout_base, Memory>::create(sub_domains, inner().getDomain(),
                                                      inner().getGhost(), inner().bc);

    // Ghost box information must be re-offloaded
    host_dev_transfer = false;
    ie_ghost<dim, domain_type, Memory, layout_base>::reset_host_dev_transfer();
  }

  grid_sm<dim, void> &getGrid() { return gr; }

  AbstractDecStrategy &inner() {
    return _inner;
  }

private:
  AbstractDecStrategy _inner;

  //! set of Boxes produced by the decomposition optimizer
  openfpm::vector<::Box<dim, size_t>> loc_box;

  //! for each sub-domain, contain the list of the neighborhood processors
  openfpm::vector<openfpm::vector<long unsigned int> > box_nn_processor;

  //! Magnification factor between distribution and
  //! decomposition
  size_t magn[dim];

  //! bool that indicate whenever the buffer has been already transfer to device
  bool host_dev_transfer = false;

  //! Processor bounding box
  ::Box<dim, domain_type> bbox;

  //! the set of all local sub-domain as vector
  SubDomains sub_domains;

  //! the remote set of all sub-domains as vector of 'sub_domains' vectors
  mutable openfpm::vector<Box_map<dim, domain_type>, Memory,
                          typename layout_base<Box_map<dim, domain_type>>::type,
                          layout_base>
      sub_domains_global;

  //! Structure that contain for each sub-sub-domain box the processor id
  //! exist for efficient global communication
  CellList<dim, domain_type, Mem_fast<Memory, int>, shift<dim, domain_type>> fine_s;

  //! Processor domain bounding box
  ::Box<dim, size_t> proc_box;

  //! Structure that decompose the space into cells without creating them
  //! useful to convert positions to CellId or sub-domain id in this case
  CellDecomposer_sm<dim, domain_type, shift<dim, domain_type>> cd;

  //! Structure that store the cartesian grid information
  grid_sm<dim, void> gr;

  //! Global sub-sub-domain graph
	Graph_CSR<nm_v<dim>, nm_e> gp;
};
#endif // SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP
