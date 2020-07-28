#ifndef OPENFPM_PDATA_MYDISTRIBUTIONSTRATEGY_UNIT_TEST_HPP
#define OPENFPM_PDATA_MYDISTRIBUTIONSTRATEGY_UNIT_TEST_HPP

#include "./utils.hpp"

#define GS_SIZE 8

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

void Parmetis_distribution_test(const unsigned int nProcs) {
  Vcluster<>& v_cl = create_vcluster();

  auto nProcUnits = v_cl.getProcessingUnits();
  if (nProcUnits != nProcs) {  // question why it was 3?
    printMe(v_cl);
    std::cerr << "# processing units = " << nProcUnits << " != 3 !"
              << std::endl;
    return;
  }

  //! [Initialize a ParMetis Cartesian graph and decompose]

  MyDistributionStrategy dist(v_cl);

  // Physical domain
  Box<3, SpaceType> box({0.0, 0.0, 0.0}, {10.0, 10.0, 10.0});

  // Grid info
  grid_sm<3, void> info({GS_SIZE, GS_SIZE, GS_SIZE});

  // Initialize Cart graph and decompose
  size_t bc[SPACE_N_DIM];
  for (size_t i = 0; i < SPACE_N_DIM; i++) {
    bc[i] = NON_PERIODIC;
  }

  dist.createCartGraph(bc, info, box);

  // First create the center of the weights distribution, check it is coherent
  // to the size of the domain
  Point<3, SpaceType> center({2.0, 2.0, 2.0});

  // It produces a sphere of radius 2.0
  // with high computation cost (5) inside the sphere and (1) outside
  setSphereComputationCosts(dist, info, center, 2.0f, 5ul, 1ul);

  // first distribution
  ParmetisGraph parmetis_graph(v_cl, v_cl.getProcessingUnits());
  dist.reset(parmetis_graph);
  parmetis_graph.decompose(dist.getVtxdist());  // decompose
  dist.distribute(parmetis_graph);              // distribute

  printMe(v_cl);
  std::cout << "assert " << parmetis_graph.get_ndec() << " == 1ul" << std::endl;

  //! [Initialize a ParMetis Cartesian graph and decompose]

  if (amIMaster(v_cl)) {
    // write the first decomposition

    dist.write("vtk_parmetis_distribution_0");
    // todo compare
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

    setSphereComputationCosts(dist, info, center, 2.0f, 5, 1);

    // With some regularity refine and write the parmetis distribution
    if ((size_t)iter % 10 == 0) {
      dist.refine(parmetis_graph);
      n_dec++;

      printMe(v_cl);
      std::cout << "assert " << parmetis_graph.get_ndec() << " == " << n_dec
                << std::endl;

      if (amIMaster(v_cl)) {
        std::stringstream str;
        str << "vtk_parmetis_distribution_" << iter;

        dist.write(str.str());
        // todo compare
      }
    }
  }

  //! [refine with parmetis the decomposition]

  printMe(v_cl);
  std::cout << "my size is " << sizeof(MyDistributionStrategy) << std::endl;
}

#endif  // OPENFPM_PDATA_MYDISTRIBUTIONSTRATEGY_UNIT_TEST_HPP
