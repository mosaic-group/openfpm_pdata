#ifndef OPENFPM_PDATA_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
#define OPENFPM_PDATA_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP

#include "./utils.hpp"

#define SUB_UNIT_FACTOR 1024

void CartDecomposition_non_periodic_test(const unsigned int nProcs) {
  Vcluster<>& vcl = create_vcluster();

  MyDecompositionStrategy dec(vcl);

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

  // Decompose
  dec.setParameters(div, box, bc, gsub);
  // breakpoint dec.decompose();
}

#endif  // OPENFPM_PDATA_MYDECOMPOSITIONSTRATEGY_UNIT_TEST_HPP
