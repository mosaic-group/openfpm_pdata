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
{
  // todo reuse AbstractDecompositionStrategy types
  using Box = SpaceBox<dim, domain_type>;

  typedef Graph_CSR<nm_v<dim>, nm_e> DecompositionGraph;

  AbstractDecStrategy _inner;

  //! Magnification factor between distribution and
  //! decomposition
  size_t magn[dim];

  //! Processor domain bounding box
  ::Box<dim, size_t> proc_box;

  //! Structure that decompose the space into cells without creating them
  //! useful to convert positions to CellId or sub-domain id in this case
  CellDecomposer_sm<dim, domain_type, shift<dim, domain_type>> cd;

  //! Structure that store the cartesian grid information
  grid_sm<dim, void> gr;

public:

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
      : _inner(v_cl) 
  {}

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

  /*! \brief Create the Cartesian graph
   */
  void createCartGraph() {
    // Create a cartesian grid graph
    CartesianGraphFactory<dim, DecompositionGraph> g_factory_part;
    inner().getGraph() = g_factory_part.template construct<NO_EDGE, nm_v_id, domain_type, dim - 1, 0>(gr.getSize(), 
                                                                                                      inner().getDomain(), 
                                                                                                      inner().bc);
  }

  void decompose() {
    createCartGraph();
  }

  void reset() {
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

  /*! \brief Return the graph of the decomposition (not-distributed)
   *
   * \return The graph
   * 
   */
  auto getGraph() -> decltype(this->_inner.getGraph())
  {
    return this->_inner.getGraph();
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

  /*! \brief Covert discrete subdomains into continuos subdomains
    *
    * \param loc_box discrete boxes
    * \param sub_domains continuos sub domains
    * 
    */
  void convertToSubDomains(openfpm::vector<::Box<dim, size_t>> & loc_box,
                           openfpm::vector<SpaceBox<dim, domain_type>,Memory,typename layout_base<SpaceBox<dim, domain_type>>::type,layout_base> & sub_domains,
                           ::Box<dim,domain_type> & bbox)
  {
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
			for (size_t i = 0 ; i < dim ; i++)
			{
				proc_box.setLow(i,0.0);
				proc_box.setHigh(i,0);

				bbox.setLow(i,0.0);
				bbox.setHigh(i,0);
			}
		}

		// convert into sub-domain
		for (size_t s = 1; s < loc_box.size(); s++)
		{
			SpaceBox<dim,domain_type> sub_d = convertDecBoxIntoSubDomain(loc_box.get(s));

			// add the sub-domain
			sub_domains.add(sub_d);

			// Calculate the bound box
			bbox.enclose(sub_d);
			proc_box.enclose(loc_box.get(s));
		}
  }


  /*! \brief Set migration cost of the vertex id
   *
   * \param id of the vertex to update
   * \param migration cost of the migration
   */
  void setMigrationCost(size_t id, size_t migration) {
    inner().getGraph().vertex(id).template get<nm_v_migration>() = migration;
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
    inner().getGraph().getChildEdge(v_id, e).template get<nm_e::communication>() =
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
    pos[0] = inner().getGraph().vertex(id).template get<nm_v_x>()[0];
    pos[1] = inner().getGraph().vertex(id).template get<nm_v_x>()[1];
    if (dim == 3) {
      pos[2] = inner().getGraph().vertex(id).template get<nm_v_x>()[2];
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
    inner().getGraph().vertex(id).template get<nm_v_computation>() = weight;
  }

  /*! \brief function that get the weight of the vertex
   * (computation cost of the sub-sub-domain id)
   *
   * \param id vertex id
   *
   */
  size_t getSubSubDomainComputationCost(size_t id) {
    return inner().getGraph().vertex(id).template get<nm_v_computation>();
  }

  /*! \brief Returns total number of sub-sub-domains in the distribution graph
   *
   * \return the total number of sub-sub-domains
   *
   */
  size_t getNSubSubDomains() { return inner().getGraph().getNVertex(); }

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
    return inner().getGraph().getNChilds(id);
  }

  grid_sm<dim, void> &getGrid() { return gr; }

  AbstractDecStrategy &inner() {
    return _inner;
  }
};
#endif // SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP
