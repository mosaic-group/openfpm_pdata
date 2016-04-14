#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Simple example
 *
 * This example show cell and verlet list of the distributed vector
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
	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();

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

		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

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
	openfpm::vector<openfpm::vector<size_t>> verlet;

	// cutting radius
	float r_cut = 1.0/(128-2);
	vd.getVerlet(verlet,r_cut);

	//
	// ### WIKI 7 ###
	//
	// for each particle we iterate across the neighborhoods particles and we
	// do some demo calculation
	//
	for (size_t i = 0 ; i < verlet.size() ; i++)
	{

		Point<3,float> p = vd.getPos(i);

		// for each neighborhood particle
		for (size_t j = 0 ; j < verlet.get(i).size() ; j++)
		{
			auto & NN = verlet.get(i);

			Point<3,float> q = vd.getPos(NN.get(j));

			// some non-sense calculation as usage demo

			// we sum the distance of all the particles
			vd.template getProp<0>(i) += p.distance(q);

			// we sum the distance of all the particles
			vd.template getProp<1>(i)[0] += p.get(0) - q.get(0);
			vd.template getProp<1>(i)[1] += p.get(0) - q.get(0);
			vd.template getProp<1>(i)[2] += p.get(0) - q.get(0);

			vd.template getProp<2>(i)[0][0] += p.get(0) - q.get(0);
			vd.template getProp<2>(i)[0][1] += p.get(0) - q.get(1);
			vd.template getProp<2>(i)[0][2] += p.get(0) - q.get(2);
			vd.template getProp<2>(i)[1][0] += p.get(1) - q.get(0);
			vd.template getProp<2>(i)[1][1] += p.get(1) - q.get(1);
			vd.template getProp<2>(i)[1][2] += p.get(1) - q.get(2);
			vd.template getProp<2>(i)[2][0] += p.get(2) - q.get(0);
			vd.template getProp<2>(i)[2][1] += p.get(2) - q.get(1);
			vd.template getProp<2>(i)[2][2] += p.get(2) - q.get(2);
		}
	}

	//
	// ### WIKI 8 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}

