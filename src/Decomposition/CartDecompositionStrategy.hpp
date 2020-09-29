#ifndef SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP

#include "Decomposition/dec_optimizer.hpp"
#include "Decomposition/AbstractDecompositionStrategy.hpp"

template <unsigned int dim, typename domain_type,
          typename AbstractDecStrategy =
              AbstractDecompositionStrategy<dim, domain_type>, typename Memory = HeapMemory, template <typename> class layout_base = memory_traits_lin>
class CartDecompositionStrategy {
  using Box = SpaceBox<dim, domain_type>;
  using DGrid = grid_sm<dim, void>;

public:
  CartDecompositionStrategy(Vcluster<> &v_cl) : dec(v_cl) {}

  ~CartDecompositionStrategy() {}

  /*! \brief Copy constructor
   *
   * \param cart object to copy
   *
   */
  CartDecompositionStrategy(
      const CartDecompositionStrategy<dim, domain_type, AbstractDecStrategy>& cart) {
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

  void setParameters(const size_t (&div_)[dim], ::Box<dim, domain_type> &domain_,
                const size_t (&bc)[dim],
                const grid_sm<dim, void> &sec_dist = grid_sm<dim, void>()) {
    std::copy(bc, bc + dim, dec.bc);

    // Set the decomposition parameters
    gr.setDimensions(div_);
    dec.domain = domain_;
    dec.cd.setDimensions(dec.domain, div_, 0);

    // calc magnification factor dec-dist
    calculate_magn(sec_dist);
  }

  /*! \brief Calculate magnification
   *
   * \param gm distribution grid
   *
   */
  void calculate_magn(const grid_sm<dim, void> &gm) {
    if (gm.size() == 0) {
      for (size_t i = 0; i < dim; i++) {
        dec.magn[i] = 1;
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

        dec.magn[i] = gr.size(i) / gm.size(i);
      }
    }
  }

  void initialize_fine_s(const ::Box<dim, domain_type> &domain) {
    dec.fine_s.clear();
    size_t div_g[dim];

    // We reduce the size of the cells by a factor 8 in 3d 4 in 2d
    for (size_t i = 0; i < dim; ++i) {
      div_g[i] = (gr.size(i) == 1) ? 1 : gr.size(i) / 2;
    }

    dec.fine_s.Initialize(domain, div_g);
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
    sub_dce.mul(dec.magn);

    // shrink by one
    for (size_t i = 0; i < dim; i++) {
      loc_box.template get<Box::p1>()[i] = sub_dce.getLow(i);
      loc_box.template get<Box::p2>()[i] = sub_dce.getHigh(i) - 1;
    }

    SpaceBox<dim, domain_type> sub_d(sub_dce);
    sub_d.mul(dec.spacing);
    sub_d += dec.domain.getP1();

    // we add the fixing sub-domains to cover all the domain
    // if (loc_box) is at the boundary we have to ensure that the box span the
    // full domain (avoiding rounding off error)
    for (size_t i = 0; i < dim; i++) {
      if (sub_dc.getHigh(i) == gr.size(i) - 1) {
        sub_d.setHigh(i, dec.domain.getHigh(i));
      }

      if (sub_dc.getLow(i) == 0) {
        sub_d.setLow(i, dec.domain.getLow(i));
      }
    }

    return sub_d;
  }

  /*! \brief Constructor, it decompose and distribute the sub-domains across the
   * processorscreateSubdomains
   *
   * \param v_cl Virtual cluster, used internally for communications
   * \param bc boundary conditions
   * \param opt option (one option is to construct)
   *
   */
  template <typename Graph>
  void createSubdomains(Graph &graph, Ghost<dim, domain_type> &ghost, DGrid gr_dist,
                        size_t opt = 0) {
    // Calculate the total number of box and and the spacing
    // on each direction
    // Get the box containing the domain
    SpaceBox<dim, domain_type> bs = dec.domain.getBox();

    for (unsigned int i = 0; i < dim; i++) {
      // Calculate the spacing
      dec.spacing[i] = (bs.getHigh(i) - bs.getLow(i)) / gr.size(i);
    }

    // fill the structure that store the processor id for each sub-domain
    initialize_fine_s(dec.domain);

    // Optimize the decomposition creating bigger spaces
    // And reducing Ghost over-stress
    dec_optimizer<dim, Graph_CSR<nm_v<dim>, nm_e>> d_o(graph,
                                                       gr_dist.getSize());

    // Ghost
    Ghost<dim, long int> ghe;

    // Set the ghost
    for (size_t i = 0; i < dim; i++)
    {
      ghe.setLow(i, static_cast<long int>(ghost.getLow(i) / dec.spacing[i]) - 1);  // -1
      ghe.setHigh(i, static_cast<long int>(ghost.getHigh(i) / dec.spacing[i]) + 1);  // +1
    }

    // optimize the decomposition or merge sub-sub-domain
    d_o.template optimize<nm_v_sub_id, nm_v_proc_id>(
        graph, dec.v_cl.getProcessUnitID(), dec.loc_box, dec.box_nn_processor, ghe, dec.bc);

    // Initialize
    if (dec.loc_box.size() > 0)
    {
      dec.bbox = convertDecBoxIntoSubDomain(dec.loc_box.get(0));
      dec.proc_box = dec.loc_box.get(0);
      dec.sub_domains.add(dec.bbox);
    }
    else
    {
      // invalidate all the boxes
      for (size_t i = 0; i < dim; i++) {
        dec.proc_box.setLow(i, 0.0);
        dec.proc_box.setHigh(i, 0);

        dec.bbox.setLow(i, 0.0);
        dec.bbox.setHigh(i, 0);
      }
    }

    // convert into sub-domain
		for (size_t s = 1; s < dec.loc_box.size(); s++)
		{
			SpaceBox<dim, domain_type> sub_d = convertDecBoxIntoSubDomain(dec.loc_box.get(s));

			// add the sub-domain
			dec.sub_domains.add(sub_d);

			// Calculate the bound box
			dec.bbox.enclose(sub_d);
			dec.proc_box.enclose(dec.loc_box.get(s));
		}

    dec.nn_prcs<dim, domain_type, layout_base, Memory>::create(dec.box_nn_processor, dec.sub_domains);
    dec.nn_prcs<dim, domain_type, layout_base, Memory>::applyBC(dec.domain, ghost, dec.bc);

    // fill fine_s structure
    // fine_s structure contain the processor id for each sub-sub-domain
    // with sub-sub-domain we mean the sub-domain decomposition before
    // running dec_optimizer (before merging sub-domains)

    construct_fine_s();  
    Initialize_geo_cell_lists();
  }

  void construct_fine_s() {
    dec.collect_all_sub_domains(dec.sub_domains_global);  // todo in this class ??

    // now draw all sub-domains in fine-s
    for (size_t i = 0; i < dec.sub_domains_global.size(); ++i) {
      // get the cells this box span
      const grid_key_dx<dim> p1 =
          dec.fine_s.getCellGrid_me(dec.sub_domains_global.template get<0>(i).getP1());
      const grid_key_dx<dim> p2 =
          dec.fine_s.getCellGrid_pe(dec.sub_domains_global.template get<0>(i).getP2());

      // Get the grid and the sub-iterator
      auto &gi = dec.fine_s.getGrid();
      grid_key_dx_iterator_sub<dim> g_sub(gi, p1, p2);

      // add the box-id to the cell list
      while (g_sub.isNext()) {
        auto key = g_sub.get();
        dec.fine_s.addCell(gi.LinId(key), i);

        ++g_sub;
      }
    }

    dec.host_dev_transfer = false;
  }

  void Initialize_geo_cell_lists() {
    // Get the processor bounding Box
    ::Box<dim, domain_type> bound = dec.getProcessorBounds();

    // Check if the box is valid
    if (bound.isValidN()) {
      // calculate the sub-divisions
      size_t div[dim];
      for (size_t i = 0; i < dim; ++i) {
        div[i] = (size_t)((bound.getHigh(i) - bound.getLow(i)) /
                          dec.cd.getCellBox().getP2()[i]);
      }

      // Initialize the geo_cell structure
      dec.ie_ghost<dim, domain_type, Memory, layout_base>::Initialize_geo_cell(bound, div);

      // Initialize shift vectors
      dec.ie_ghost<dim, domain_type, Memory, layout_base>::generateShiftVectors(dec.domain, dec.bc);
    }
  }

  template <typename Model, typename Graph>
  void decompose(Model m, Graph &graph, openfpm::vector<rid> &vtxdist) {
    dec.decompose(m, graph, vtxdist);
  }

  template <typename Graph>
  void merge(Graph &graph, Ghost<dim, domain_type> &ghost, DGrid gr_dist) {
    createSubdomains(graph, ghost, gr_dist);
    dec.calculateGhostBoxes(ghost);
  }

  DGrid &getGrid() { return gr; }

// todo private:
  AbstractDecStrategy dec;
  
  //! Structure that store the cartesian grid information
  DGrid gr;
};
#endif // SRC_DECOMPOSITION_CART_DECOMPOSITION_STRATEGY_HPP
