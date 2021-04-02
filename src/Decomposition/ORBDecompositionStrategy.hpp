#ifndef SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP
#define SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP

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
#include "Decomposition/ORB.hpp"

template <unsigned int dim, typename domain_type, typename Memory = HeapMemory, template <typename> class layout_base = memory_traits_lin,
          typename AbstractDecStrategy =
              AbstractDecompositionStrategy<dim, domain_type>>
class OrbDecompositionStrategy {

  using Box = SpaceBox<dim, domain_type>;
  using Orb = ORB<dim, domain_type>;

  typedef openfpm::vector<SpaceBox<dim, domain_type>, Memory, typename layout_base<SpaceBox<dim, domain_type>>::type,layout_base> vector_subdomains_type;

  void go_deep(size_t v_id, ::Box<dim,domain_type> box, typename Orb::graph_type & graph, vector_subdomains_type & sub_domains)
  {
    // Ho childe
    if (graph.getNChilds(v_id) != 0)  // todo 08.04 fix 2 boxes
    {
      for (int i = 0 ; i < graph.getNChilds(v_id) ; i++)
      {
        ::Box<dim,domain_type> bc = box;
        int dir = graph.template vertex_p<ORB_node<domain_type>::dir_split>(v_id);
        domain_type CM = graph.template vertex_p<ORB_node<domain_type>::CM>(v_id);

          if (i == 0)
          {box.setHigh(dir,CM);}

          if (i == 1)
          {box.setLow(dir,CM);}
          go_deep(graph.getChild(v_id,i),box,graph,sub_domains);
      }
    } else {  // I'm leaf
      sub_domains.add(box);
    }
  }

public:
  OrbDecompositionStrategy(Vcluster<> &v_cl) : _inner(v_cl) {}

  /*! \brief Copy constructor
   *
   * \param cart object to copy
   *
   */
  OrbDecompositionStrategy(
      const OrbDecompositionStrategy<dim, domain_type, AbstractDecStrategy>& cart) {
    this->operator=(cart);
  }

  /*! \brief Copy the element
   *
   * \param cart element to copy
   *
   * \return itself
   *
   */
  OrbDecompositionStrategy<dim, domain_type> &
  operator=(const AbstractDecStrategy &cart) {
    // todo
  }

  ~OrbDecompositionStrategy() {
    if (orb) {
      delete orb;
    }
  }


  void setParameters(const size_t (&div_)[dim],
    ::Box<dim, domain_type> &domain_,
    const size_t (&bc)[dim],
    const Ghost<dim, domain_type> &ghost,
    const grid_sm<dim, void> &sec_dist = grid_sm<dim, void>()) {
    inner().setDomain(domain_);
  }

  /*! \brief Stub method to homogenize the interface do not use
   *
   */
  std::pair<float, size_t> computeCommunicationCosts() {
  }

  /*! \brief Stub method to homogenize the interface do not use
   *
   */
  void createCartGraph() {
  }

  void decompose(openfpm::vector<Point<dim, domain_type>> points) {
    orb = new Orb(inner().getDomain(), inner().getVcluster().getProcessingUnits(), points);
  }

    /*! \brief Covert discrete subdomains into continuos subdomains
    *
    * \param loc_box discrete boxes
    * \param sub_domains continuos sub domains
    * 
    */
  void convertToSubDomains(openfpm::vector<::Box<dim, size_t>> & loc_box,
                           vector_subdomains_type & sub_domains,
                           ::Box<dim,domain_type> & bbox) {
    // la bbox e' la scatola che contiene tutte le scatole

    auto & graph = orb->getGraph();

    go_deep(0,inner().getDomain(),graph,sub_domains);

    // enclose all the rest -> convert into sub-domain
    for (auto s = 1; s < loc_box.size(); ++s) {
        //SpaceBox<dim, domain_type> sub_d = loc_box.get(s);
        //sub_domains.add(sub_d);  // add the sub-domain
        //bbox.enclose(sub_d);  // Calculate the bound box
    }
  }

  Graph_CSR<nm_v<dim>, nm_e> &getGraph() {
    //orb -> graph (V box)


    // todo get leaves from auto tree = orb->grp;  // Graph_CSR<ORB_node<T>, no_edge>
    // todo test only: ~trivial decomposition

    // build it
    if (inner().getGraph().getNVertex() < inner().v_cl.getProcessingUnits()) {
      // create only necessary vertices in the graph
      for (auto i = inner().getGraph().getNVertex(); i < inner().v_cl.getProcessingUnits(); ++i) {
        inner().getGraph().addVertex();
      }
    }

    return inner().getGraph();
  }

  AbstractDecStrategy &inner() {
    return _inner;
  }

private:
  AbstractDecStrategy _inner;
  Orb* orb;
};
#endif // SRC_DECOMPOSITION_ORB_DECOMPOSITION_STRATEGY_HPP
