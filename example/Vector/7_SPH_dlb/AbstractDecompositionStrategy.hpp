#ifndef OPENFPM_PDATA_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
#define OPENFPM_PDATA_ABSTRACT_DECOMPOSITION_STRATEGY_HPP

#include "config.h"
#include <cmath>
#include "VCluster/VCluster.hpp"
#include "Graph/CartesianGraphFactory.hpp"
#include "Decomposition/Decomposition.hpp"
#include "Vector/map_vector.hpp"
#include <vector>
#include <initializer_list>
#include "SubdomainGraphNodes.hpp"
#include "Decomposition/dec_optimizer.hpp"
#include "Space/Shape/Box.hpp"
#include "Space/Shape/Point.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include <unordered_map>
#include "NN/CellList/CellList.hpp"
#include "Space/Ghost.hpp"
#include "Decomposition/common.hpp"
#include "Decomposition/ie_loc_ghost.hpp"
#include "Decomposition/ie_ghost.hpp"
#include "Decomposition/nn_processor.hpp"
#include "GraphMLWriter/GraphMLWriter.hpp"
#include "DLB/DLB.hpp"
#include "util/se_util.hpp"
#include "util/mathutil.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/Domain_NN_calculator_cart.hpp"
#include "Decomposition/Domain_icells_cart.hpp"

template<unsigned int dim, typename T, typename Memory = HeapMemory, template <typename> class layout_base = memory_traits_lin>
class AbstractDecompositionStrategy: public ie_loc_ghost<dim, T,layout_base, Memory>,
                         public nn_prcs<dim, T,layout_base,Memory>,
                         public ie_ghost<dim,T,Memory,layout_base>,
                         public domain_nn_calculator_cart<dim>,
                         public domain_icell_calculator<dim,T,layout_base,Memory> {
public:

    //! Runtime virtual cluster machine
    Vcluster<> & v_cl;  // question can be private?

    //! Type of the domain we are going to decompose
    using domain_type = T;

    //! It simplify to access the SpaceBox element
    using Box = SpaceBox<dim, T>;

    /*! \brief Abstract decomposition constructor
     *
     * \param v_cl Virtual cluster, used internally to handle or pipeline communication
     *
     */
    AbstractDecompositionStrategy(Vcluster<> & v_cl) :nn_prcs<dim, T, layout_base,Memory>(v_cl), v_cl(v_cl) {}

    /*! \brief Function that set the computational cost for a of a sub-sub domain
	 *
	 * \param id vertex id
	 * \param weight compotational cost
	 *
	 */
    inline void setSubSubDomainComputationCost(size_t id, size_t weight) {}
    
    /*! \brief Add computation cost i to the subsubdomain with global id gid
	 *
	 * \param gid global id of the subsubdomain to update
	 * \param i Cost increment
	 */
    inline void addComputationCost(size_t gid, size_t i) {}

    /*! \brief function that return the computation cost of the sub-sub-domain id
	 *
	 * \param id sub-sub-domain id
	 *
	 * \return the computational cost
	 *
	 */
    inline size_t getSubSubDomainComputationCost(size_t id) {
        return 0;
    }

    /*! \brief Calculate communication and migration costs
	 *
	 * \param ts how many timesteps have passed since last calculation, used to approximate the cost
	 */
    void computeCommunicationAndMigrationCosts(size_t ts) {}
    
    /*! \brief Return the box of the physical domain
	 *
	 * \return The physical domain box
	 *
	 */
    const ::Box<dim,T> & getDomain() const
    {
      return domain;
    }
    
    /*! \brief Distribution grid
	 *
	 * \return the grid
	 */
    const grid_sm<dim, void> getDistGrid()
    {
      return gr_dist;
    }

    /*! \brief Start decomposition
	 *
	 */
    template<typename Model>
    void decompose(Model m) {}

protected:
    //! rectangular domain to decompose
    ::Box<dim,T> domain;

    //! Structure that store the cartesian grid information
    grid_sm<dim, void> gr_dist;
};
#endif //OPENFPM_PDATA_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
