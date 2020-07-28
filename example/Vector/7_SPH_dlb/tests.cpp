#include <iostream>

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"

#include "SubdomainGraphNodes.hpp"
#include "Decomposition/Distribution/parmetis_util.hpp"
#include "Graph/ids.hpp"
#include "Graph/CartesianGraphFactory.hpp"

#include "./MyStuff.hpp"

#define GS_SIZE 8

bool amIMaster(Vcluster<>& v_cl) { return v_cl.getProcessUnitID() == 0; }

template <unsigned int dim, typename Distribution>
void setSphereComputationCosts(Distribution& dist,
                               grid_sm<dim, void>& gr,
                               Point<3, SpaceType> center,
                               SpaceType radius,
                               size_t max_l,
                               size_t min_l) {
  SpaceType radius2 = radius * radius;
  SpaceType eq;

  // Position structure for the single vertex
  SpaceType pos[dim];

  for (size_t i = 0; i < dist.getNSubSubDomains(); i++) {
    dist.getSubSubDomainPosition(i, pos);

    eq = 0;
    for (size_t j = 0; j < dim; j++)
      eq += (pos[j] - center.get(j)) * (pos[j] - center.get(j));

    if (eq <= radius2) {
      dist.setComputationCost(i, max_l);
      dist.setMigrationCost(i, max_l * 2);
    } else {
      dist.setComputationCost(i, min_l);
      dist.setMigrationCost(i, min_l * 2);
    }

    // set Migration cost and communication cost
    for (size_t j = 0; j < dist.getNSubSubDomainNeighbors(i); j++) {
      dist.setCommunicationCost(i, j, 1);
    }
  }
}

void doTest() {
  Vcluster<>& v_cl = create_vcluster();

  auto nProcUnits = v_cl.getProcessingUnits();
  if (nProcUnits != 3) {  // question why 3?
    std::cerr << "# processing units = " << nProcUnits << " != 3 !"
              << std::endl;
    return;
  }

  //! [Initialize a ParMetis Cartesian graph and decompose]
  MyDistributionStrategy pmet_dist(v_cl);

  // Physical domain
  Box<3, SpaceType> box({0.0, 0.0, 0.0}, {10.0, 10.0, 10.0});

  // Grid info
  grid_sm<3, void> info({GS_SIZE, GS_SIZE, GS_SIZE});

  // Initialize Cart graph and decompose
  // todo pmet_dist.createCartGraph(info, box);

  // First create the center of the weights distribution, check it is coherent
  // to the size of the domain
  Point<3, SpaceType> center({2.0, 2.0, 2.0});

  // It produces a sphere of radius 2.0
  // with high computation cost (5) inside the sphere and (1) outside
  setSphereComputationCosts(pmet_dist, info, center, 2.0f, 5ul, 1ul);

  // first decomposition
  // todo pmet_dist.decompose();

  // todo BOOST_REQUIRE_EQUAL(pmet_dist.get_ndec(), 1ul);

  //! [Initialize a ParMetis Cartesian graph and decompose]

  if (amIMaster(v_cl)) {
    // write the first decomposition
    // todo pmet_dist.write("vtk_parmetis_distribution_0");
    // todo compare
    // todo BOOST_REQUIRE_EQUAL(true, test);
  }

  //! [refine with parmetis the decomposition]

  SpaceType stime = 0.0, etime = 10.0, tstep = 0.1;

  // Shift of the sphere at each iteration
  Point<3, SpaceType> shift({tstep, tstep, tstep});

  size_t iter = 1;
  size_t n_dec = 1;

  for (SpaceType t = stime; t < etime; t = t + tstep, iter++) {
    if (t < etime / 2) {
      center += shift;
    } else {
      center -= shift;
    }

    setSphereComputationCosts(pmet_dist, info, center, 2.0f, 5, 1);

    // With some regularity refine and write the parmetis distribution
    if ((size_t)iter % 10 == 0) {
      // todo pmet_dist.refine();
      n_dec++;
      // todo BOOST_REQUIRE_EQUAL(pmet_dist.get_ndec(), n_dec);

      if (amIMaster(v_cl)) {
        std::stringstream str;
        str << "vtk_parmetis_distribution_" << iter;

        // todo pmet_dist.write(str.str());
        // todo compare
        // todo BOOST_REQUIRE_EQUAL(true, test);
      }
    }
  }

  //! [refine with parmetis the decomposition]

  // todo BOOST_REQUIRE_EQUAL(sizeof(ParMetisDistribution<3,SpaceType>),872ul);
}

int main(int argc, char* argv[]) {
  openfpm_init(&argc, &argv);

  doTest();

  openfpm_finalize();
}