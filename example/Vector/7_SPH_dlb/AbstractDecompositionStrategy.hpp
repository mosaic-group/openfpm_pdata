#ifndef OPENFPM_PDATA_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
#define OPENFPM_PDATA_ABSTRACT_DECOMPOSITION_STRATEGY_HPP

#include <cmath>
#include <initializer_list>
#include <unordered_map>
#include <utility>
#include <vector>
#include "DLB/DLB.hpp"
#include "Decomposition/Decomposition.hpp"
#include "Decomposition/Domain_icells_cart.hpp"
#include "Decomposition/common.hpp"
#include "Decomposition/dec_optimizer.hpp"
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

template <unsigned int dim,
          typename T,
          typename Memory = HeapMemory,
          template <typename> class layout_base = memory_traits_lin>
class AbstractDecompositionStrategy
  : public ie_loc_ghost<dim, T, layout_base, Memory>,
    public nn_prcs<dim, T, layout_base, Memory>,
    public ie_ghost<dim, T, Memory, layout_base>,
    public domain_icell_calculator<dim, T, layout_base, Memory> {
  //! Type of the domain we are going to decompose
  using domain_type = T;

  //! It simplify to access the SpaceBox element
  using Box = SpaceBox<dim, T>;

public:
  //! Structure that decompose the space into cells without creating them
  //! useful to convert positions to CellId or sub-domain id in this case
  CellDecomposer_sm<dim, T, shift<dim, T>> cd;  // todo private

  //! ghost info
  Ghost<dim, T> ghost;  // todo private

  /*! \brief Abstract decomposition constructor
   *
   * \param v_cl Virtual cluster, used internally to handle or pipeline
   * communication
   *
   */
  AbstractDecompositionStrategy(Vcluster<>& v_cl)
    : nn_prcs<dim, T, layout_base, Memory>(v_cl), v_cl(v_cl) {}

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
  inline size_t getSubSubDomainComputationCost(size_t id) { return 0; }

  /*! \brief Calculate communication and migration costs
   *
   * \param ts how many timesteps have passed since last calculation, used to
   * approximate the cost
   */
  void computeCommunicationAndMigrationCosts(size_t ts) {}

  std::pair<float, size_t> computeCommunicationCosts() {
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

    return std::make_pair(migration, norm);
  }

  /*! \brief Return the box of the physical domain
   *
   * \return The physical domain box
   *
   */
  const ::Box<dim, T>& getDomain() const { return domain; }

  /*! \brief Distribution grid
   *
   * \return the grid
   */
  const grid_sm<dim, void> getDistGrid() { return gr_dist; }

  unsigned int getDim() const { return dim; }

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

  /*! \brief Return the ghost
   *
   *
   * \return the ghost extension
   *
   */
  const Ghost<dim, T>& getGhost() const { return ghost; }

  /*! \brief Start decomposition
   *
   */
  template <typename Model>
  void decompose(Model m) {}

  void merge() {
    createSubdomains();
    calculateGhostBoxes();
  }

  void onEnd() {
    domain_icell_calculator<dim, T, layout_base, Memory>::
        CalculateInternalCells(
            v_cl,
            ie_ghost<dim, T, Memory, layout_base>::private_get_vb_int_box(),
            sub_domains,
            this->getProcessorBounds(),
            this->getGhost().getRcut(),
            this->getGhost());
  }

protected:
  //! rectangular domain to decompose
  ::Box<dim, T> domain;

  //! Processor bounding box
  ::Box<dim, T> bbox;

  //! Structure that store the cartesian grid information
  grid_sm<dim, void> gr_dist;

private:
  //! Runtime virtual cluster machine
  Vcluster<>& v_cl;  // question can be private?

  //! the set of all local sub-domain as vector
  openfpm::vector<Box, Memory, typename layout_base<Box>::type, layout_base>
      sub_domains;

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

  /*! \brief Constructor, it decompose and distribute the sub-domains across the
   * processors
   *
   * \param v_cl Virtual cluster, used internally for communications
   * \param bc boundary conditions
   * \param opt option (one option is to construct)
   *
   */
  void createSubdomains(size_t opt = 0) {}

  /*! \brief It calculate the internal ghost boxes
   *
   * Example: Processor 10 calculate
   * B8_0 B9_0 B9_1 and B5_0
   *
   *
   *
   \verbatim

  +----------------------------------------------------+
  |                                                    |
  |                 Processor 8                        |
  |                 Sub+domain 0 +-----------------------------------+ | | | |
|                                   |
  ++--------------+---+---------------------------+----+        Processor 9 | |
|   |     B8_0                  |    |        Subdomain 0                | |
+------------------------------------+                                   | | |
|                           |    |                                   | | |   |
|B9_0|                                   | |              | B |    Local
processor        |    |                                   | | Processor 5  | 5 |
Subdomain 0            |    |                                   | | Subdomain 0
| _ |                           +----------------------------------------+ | | 0
|                           |    |                                   | | |   |
|    |                                   | |              |   | |    | Processor
9                | |              |   |                           |B9_1|
Subdomain 1                | |              |   |                           | |
| |              |   |                           |    | | |              |   |
|    |                                   |
   +--------------+---+---------------------------+----+ | | |
                             +-----------------------------------+


 \endverbatim

       and also
       G8_0 G9_0 G9_1 G5_0 (External ghost boxes)

\verbatim

      +----------------------------------------------------+
      |                 Processor 8                        |
      |                 Subdomain 0 +-----------------------------------+ | | |
      |           +---------------------------------------------+ | | | G8_0 |
|                              |
  +-----+---------------+------------------------------------+    |   Processor
9                | |                 |   |                                    |
|   Subdomain 0                | |                 |   | |G9_0| | | |   | |    |
| |                 |   |                                    |    | | | |   |
Local processor             |    |                              | |  Processor 5
|   |        Sub+domain 0                |    |                              |
  |  Subdomain 0    |   | +-----------------------------------+ | |   | |    | |
  |                 | G |                                    |    | | | | 5 | |
|   Processor 9                | |                 | | | |    |   Subdomain 1 |
  |                 | 0 |                                    |G9_1| | | |   | |
|                              | |                 |   | |    | |
  +---------------------+------------------------------------+    | | | |    | |
            +----------------------------------------+----+------------------------------+

   \endverbatim

   *
   * ghost margins for each dimensions (p1 negative part) (p2 positive part)
   *
   *
   \verbatim

               ^ p2[1]
               |
               |
           +----+----+
           |         |
           |         |
   p1[0]<-----+         +----> p2[0]
           |         |
           |         |
           +----+----+
               |
               v  p1[1]

   \endverbatim

   *
   *
   */
  void calculateGhostBoxes() {}
};
#endif  // OPENFPM_PDATA_ABSTRACT_DECOMPOSITION_STRATEGY_HPP
