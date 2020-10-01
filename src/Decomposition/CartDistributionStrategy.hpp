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

    dist.distribute(graph);
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
