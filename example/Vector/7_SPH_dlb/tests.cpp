#include <iostream>

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"

#include "SubdomainGraphNodes.hpp"
#include "Decomposition/Distribution/parmetis_util.hpp"
#include "Graph/ids.hpp"
#include "Graph/CartesianGraphFactory.hpp"

#include "./MyStuff.hpp"
#include "./MyDecompositionStrategy_unit_test.hpp"

int main(int argc, char* argv[]) {
  openfpm_init(&argc, &argv);

  CartDecomposition_non_periodic_test(2);

  openfpm_finalize();
}