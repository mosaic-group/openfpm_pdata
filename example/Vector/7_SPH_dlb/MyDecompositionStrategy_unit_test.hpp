#ifndef OPENFPM_PDATA_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
#define OPENFPM_PDATA_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP

#include "./utils.hpp"

#define SUB_UNIT_FACTOR 1024

void CartDecomposition_non_periodic_test(const unsigned int nProcs) {
  Vcluster<>& vcl = create_vcluster();

  // specify
  // - how we want to add the computational cost ...
  MyComputationalCostsModel mcc;

  // - how we want to decompose ...
  MyDecompositionStrategy dec(vcl);
  MyDecompositionModel mde;

  // - how we want to distribute ...
  MyDistributionStrategy dist(vcl);
  MyDistributionModel mdi;

  // ... and our shared information

  //! Convert the graph to parmetis format
  ParmetisGraph parmetis_graph(vcl, vcl.getProcessingUnits());

  // Physical domain
  Box<3, SpaceType> box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
  size_t div[3];
  size_t div_sub[3];

  // Get the number of processor and calculate the number of sub-domain
  // for each processor (SUB_UNIT_FACTOR=64)
  size_t n_proc = vcl.getProcessingUnits();
  size_t n_sub = n_proc * SUB_UNIT_FACTOR * 4 * 4 * 4;

  // Set the number of sub-domains on each dimension (in a scalable way)
  for (int i = 0; i < 3; i++) {
    div[i] = openfpm::math::round_big_2(pow(n_sub, 1.0 / 3));
  }

  // create a sub_distribution grid
  for (int i = 0; i < 3; i++) {
    div_sub[i] = div[i] / 4;
  }

  grid_sm<3, void> gsub(div_sub);

  // Define ghost
  Ghost<3, SpaceType> g(0.01);

  // Boundary conditions
  size_t bc[] = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

  // todo why don't we compute comm costs?

  // init
  dec.setParameters(div, box, bc, gsub);
  dist.setParameters(dec.getGrid(), gsub);
  dist.createCartGraph(bc, box);

  //////////////////////////////////////////////////////////////////// decompose
  dec.reset();
  dist.reset(parmetis_graph);
  if (dec.shouldSetCosts()) {
    mcc.computeCommunicationAndMigrationCosts(dec, dist);
  }
  dec.decompose(mde, parmetis_graph, dist.getVtxdist());

  /////////////////////////////////////////////////////////////////// distribute
  dist.distribute(parmetis_graph);

  ///////////////////////////////////////////////////////////////////////  merge
  dec.merge(dist.getGraph(), dist.getGhost(), dist.getGrid());

  ///////////////////////////////////////////////////////////////////// finalize
  dist.onEnd(dec.getSubDomains());
  dec.onEnd(dist.getGhost());

  // For each calculated ghost box
  printVar(dec.getNIGhostBox());
  for (size_t i = 0; i < dec.getNIGhostBox(); i++) {
    SpaceBox<3, float> b = dec.getIGhostBox(i);
    size_t proc = dec.getIGhostBoxProcessor(i);

    // sample one point inside the box
    Point<3, float> p = b.rnd();

    // Check that ghost_processorsID return that processor number
    const openfpm::vector<size_t>& pr =
        dec.ghost_processorID<CartDecomposition<3, float>::processor_id>(p);

    bool found = isIn(pr, proc);
    if (!found) {
      const openfpm::vector<size_t> pr2 =
          dec.ghost_processorID<CartDecomposition<3, float>::processor_id>(p);
    }

    printMe(vcl);
    std::cout << "assert " << found << " == true" << std::endl;
  }

  // Check the consistency
  bool val = dec.check_consistency();

  printMe(vcl);
  std::cout << "assert " << val << " == true" << std::endl;
}

#endif  // OPENFPM_PDATA_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
