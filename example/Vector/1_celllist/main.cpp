
#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 *
 * This example show cell lists for the distributed vector
 *
 * ### WIKI END ###
 *
 */

int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library, we create a Box that define our domain, boundary conditions, ghost
	// and the grid size
	//
	init_global_v_cluster(&argc,&argv);
	Vcluster & v_cl = *global_v_cluster;

	// we create a 128x128x128 Grid iterator
	size_t sz[3] = {128,128,128};

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(1.0/(128-2));

	//
	// ### WIKI 3 ###
	//
	// Here we define a distributed vector in 3D, containing 3 properties, a
	// scalar double, a vector double[3], and a tensor or rank 2 double[3][3].
	// In this case the vector contain 0 particles in total
	//
	vector_dist<3,float, aggregate<double,double[3],double[3][3]> > vd(0,box,bc,ghost);

	//
	// ### WIKI 4 ###
	//
	// We define a grid iterator, to create particles on a grid like way.
	// An important note is that the grid iterator, iterator only on the
	// local nodes for each processor for example suppose to have a domain like
	// the one in figure
	//
	//   +---------+
	//   |* * *|* *|
	//   |  2  |   |
	//   |* * *|* *|
	//   |   ---   |
	//   |* *|* * *|
	//   |   |     |
	//   |* *|* * *|
	//   |   |  1  |
	//   |* *|* * *|
	//   +---------+
	//
	// divided in 2 processors, the processor 1 will iterate only on the points
	// inside the portion of space marked with one. A note grid iterator follow the
	// boundary condition specified in vector. For a perdiodic 2D 5x5 grid we have
	//
	//   +---------+
	//   * * * * * |
	//   |         |
	//   * * * * * |
	//   |         |
	//   * * * * * |
	//   |         |
	//   * * * * * |
	//   |         |
	//   *-*-*-*-*-+
	//
	// Because the right border is equivalent to the left border, while for a non periodic we have the
	// following distribution of points
	//
	//   *-*-*-*-*
	//   |       |
	//   * * * * *
	//   |       |
	//   * * * * *
	//   |       |
	//   * * * * *
	//   |       |
	//   *-*-*-*-*
	//
	// So in this loop each processor will place particles on a grid
	//
	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		vd.add();

		auto key = it.get();

		vd.template getLastPos<0>()[0] = key.get(0) * it.getSpacing(0);
		vd.template getLastPos<0>()[1] = key.get(1) * it.getSpacing(1);
		vd.template getLastPos<0>()[2] = key.get(2) * it.getSpacing(2);

		++it;
	}

	//
	// ### WIKI 5 ###
	//
	// we synchronize the ghost, the scalar property, the vector, and the rank 2 tensor
	// (just for fun)
	vd.ghost_get<0,1,2>();

	//
	// ### WIKI 6 ###
	//
	// If the particle does not move, or does not move that much we can create a verlet list
	// for each particle, it internally use CellList to find the neighborhood but it is still
	// an expensive operation
	//
	float r_cut = 1.0/(128-2);
	auto NN = vd.getCellList(r_cut);

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		Point<3,float> xp = vd.getPos<0>(p);

		auto Np = NN.getIterator(NN.getCell(vd.getPos<0>(p)));

		while (Np.isNext())
		{
			auto q = Np.get();

			// repulsive

			Point<3,float> xq = vd.getPos<0>(q);
			Point<3,float> f = (xp - xq);

			// we sum the distance of all the particles
			vd.template getProp<0>(p) += f.norm();;

			// we sum the distance of all the particles
			vd.template getProp<1>(p)[0] += f.get(0);
			vd.template getProp<1>(p)[1] += f.get(0);
			vd.template getProp<1>(p.getKey())[2] += f.get(0);

			vd.template getProp<2>(p)[0][0] += xp.get(0) - xq.get(0);
			vd.template getProp<2>(p)[0][1] += xp.get(0) - xq.get(1);
			vd.template getProp<2>(p)[0][2] += xp.get(0) - xq.get(2);
			vd.template getProp<2>(p)[1][0] += xp.get(1) - xq.get(0);
			vd.template getProp<2>(p)[1][1] += xp.get(1) - xq.get(1);
			vd.template getProp<2>(p)[1][2] += xp.get(1) - xq.get(2);
			vd.template getProp<2>(p)[2][0] += xp.get(2) - xq.get(0);
			vd.template getProp<2>(p)[2][1] += xp.get(2) - xq.get(1);
			vd.template getProp<2>(p)[2][2] += xp.get(2) - xq.get(2);
		}

		++it2;
	}

	//
	// ### WIKI 10 ###
	//
	// Deinitialize the library
	//
	delete_global_v_cluster();
}




