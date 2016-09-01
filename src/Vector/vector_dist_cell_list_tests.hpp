/*
 * vector_dist_cell_list_tests.hpp
 *
 *  Created on: Aug 16, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_CELL_LIST_TESTS_HPP_
#define SRC_VECTOR_VECTOR_DIST_CELL_LIST_TESTS_HPP_


///////////////////////// test hilb ///////////////////////////////

BOOST_AUTO_TEST_CASE( vector_dist_reorder_2d_test )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 48)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 524288 * v_cl.getProcessingUnits();

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v( "Testing 2D vector with hilbert curve reordering k<=",k);

	// 2D test
	for ( ; k >= 2 ; k-= decrement(k,big_step) )
	{
		BOOST_TEST_CHECKPOINT( "Testing 2D vector with hilbert curve reordering k=" << k );

		Box<2,float> box({0.0,0.0},{1.0,1.0});

		// Boundary conditions
		size_t bc[2]={NON_PERIODIC,NON_PERIODIC};

		vector_dist<2,float, Point_test<float>, CartDecomposition<2,float> > vd(k,box,bc,Ghost<2,float>(0.0));

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);

			++it;
		}

		vd.map();

		// Create first cell list

		auto NN1 = vd.getCellList(0.01);

		//An order of a curve
		int32_t m = 6;

		//Reorder a vector
		vd.reorder(m);

		// Create second cell list
		auto NN2 = vd.getCellList(0.01);

		//Check equality of cell sizes
		for (size_t i = 0 ; i < NN1.getGrid().size() ; i++)
		{
			size_t n1 = NN1.getNelements(i);
			size_t n2 = NN2.getNelements(i);

			BOOST_REQUIRE_EQUAL(n1,n2);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_cl_random_vs_hilb_forces_test )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 48)
		return;

	///////////////////// INPUT DATA //////////////////////

	// Dimensionality of the space
	const size_t dim = 3;
	// Cut-off radiuses. Can be put different number of values
	openfpm::vector<float> cl_r_cutoff {0.05};
	// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
	size_t cl_k_start = 10000;
	// The lower threshold for number of particles
	size_t cl_k_min = 1000;
	// Ghost part of distributed vector
	double ghost_part = 0.05;

	///////////////////////////////////////////////////////

	//For different r_cut
	for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
	{
		//Cut-off radius
		float r_cut = cl_r_cutoff.get(r);

		//Number of particles
		size_t k = cl_k_start * v_cl.getProcessingUnits();

		std::string str("Testing " + std::to_string(dim) + "D vector's forces (random vs hilb celllist) k<=");

		vector_dist_test::print_test_v(str,k);

		//For different number of particles
		for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
		{
			BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector's forces (random vs hilb celllist) k<=" << k_int );

			Box<dim,float> box;

			for (size_t i = 0; i < dim; i++)
			{
				box.setLow(i,0.0);
				box.setHigh(i,1.0);
			}

			// Boundary conditions
			size_t bc[dim];

			for (size_t i = 0; i < dim; i++)
				bc[i] = PERIODIC;

			vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(ghost_part));

			vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd2(k_int,box,bc,Ghost<dim,float>(ghost_part));

			// Initialize dist vectors
			vd_initialize_double<dim>(vd, vd2, v_cl, k_int);

			vd.template ghost_get<0>();
			vd2.template ghost_get<0>();

			//Get a cell list

			auto NN = vd.getCellList(r_cut);

			//Calculate forces

			calc_forces<dim>(NN,vd,r_cut);

			//Get a cell list hilb

			auto NN_hilb = vd2.getCellList_hilb(r_cut);

			//Calculate forces
			calc_forces_hilb<dim>(NN_hilb,vd2,r_cut);

			// Calculate average
			size_t count = 1;
			Point<dim,float> avg;
			for (size_t i = 0 ; i < dim ; i++)	{avg.get(i) = 0.0;}

			auto it_v2 = vd.getIterator();
			while (it_v2.isNext())
			{
				//key
				vect_dist_key_dx key = it_v2.get();

				for (size_t i = 0; i < dim; i++)
					avg.get(i) += fabs(vd.template getProp<0>(key)[i]);

				++count;
				++it_v2;
			}

			for (size_t i = 0 ; i < dim ; i++)	{avg.get(i) /= count;}

			auto it_v = vd.getIterator();
			while (it_v.isNext())
			{
				//key
				vect_dist_key_dx key = it_v.get();

				for (size_t i = 0; i < dim; i++)
				{
					auto a1 = vd.template getProp<0>(key)[i];
					auto a2 = vd2.template getProp<0>(key)[i];

					//Check that the forces are (almost) equal
					float per = 0.1;
					if (a1 != 0.0)
						per = fabs(0.1*avg.get(i)/a1);

					BOOST_REQUIRE_CLOSE(a1,a2,per);
				}

				++it_v;
			}
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_cl_random_vs_reorder_forces_test )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 48)
		return;

	///////////////////// INPUT DATA //////////////////////

	// Dimensionality of the space
	const size_t dim = 3;
	// Cut-off radiuses. Can be put different number of values
	openfpm::vector<float> cl_r_cutoff {0.01};
	// The starting amount of particles (remember that this number is multiplied by number of processors you use for testing)
	size_t cl_k_start = 10000;
	// The lower threshold for number of particles
	size_t cl_k_min = 1000;
	// Ghost part of distributed vector
	double ghost_part = 0.01;

	///////////////////////////////////////////////////////

	//For different r_cut
	for (size_t r = 0; r < cl_r_cutoff.size(); r++ )
	{
		//Cut-off radius
		float r_cut = cl_r_cutoff.get(r);

		//Number of particles
		size_t k = cl_k_start * v_cl.getProcessingUnits();

		std::string str("Testing " + std::to_string(dim) + "D vector's forces (random vs reorder) k<=");

		vector_dist_test::print_test_v(str,k);

		//For different number of particles
		for (size_t k_int = k ; k_int >= cl_k_min ; k_int/=2 )
		{
			BOOST_TEST_CHECKPOINT( "Testing " << dim << "D vector's forces (random vs reorder) k<=" << k_int );

			Box<dim,float> box;

			for (size_t i = 0; i < dim; i++)
			{
				box.setLow(i,0.0);
				box.setHigh(i,1.0);
			}

			// Boundary conditions
			size_t bc[dim];

			for (size_t i = 0; i < dim; i++)
				bc[i] = PERIODIC;

			vector_dist<dim,float, aggregate<float[dim], float[dim]>, CartDecomposition<dim,float> > vd(k_int,box,bc,Ghost<dim,float>(ghost_part));

			// Initialize vd
			vd_initialize<dim,decltype(vd)>(vd, v_cl, k_int);

			vd.template ghost_get<0>();

			//Get a cell list

			auto NN1 = vd.getCellList(r_cut);

			//Calculate forces '0'

			calc_forces<dim>(NN1,vd,r_cut);

			//Reorder and get a cell list again

			vd.reorder(4);

			vd.template ghost_get<0>();

			auto NN2 = vd.getCellList(r_cut);

			//Calculate forces '1'
			calc_forces<dim,1>(NN2,vd,r_cut);

			// Calculate average (For Coverty scan we start from 1)
			size_t count = 1;
			Point<dim,float> avg;
			for (size_t i = 0 ; i < dim ; i++)	{avg.get(i) = 0.0;}

			auto it_v2 = vd.getIterator();
			while (it_v2.isNext())
			{
				//key
				vect_dist_key_dx key = it_v2.get();

				for (size_t i = 0; i < dim; i++)
					avg.get(i) += fabs(vd.template getProp<0>(key)[i]);

				++count;
				++it_v2;
			}

			for (size_t i = 0 ; i < dim ; i++)	{avg.get(i) /= count;}

			//Test for equality of forces
			auto it_v = vd.getDomainIterator();

			while (it_v.isNext())
			{
				//key
				vect_dist_key_dx key = it_v.get();

				for (size_t i = 0; i < dim; i++)
				{
					auto a1 = vd.template getProp<0>(key)[i];
					auto a2 = vd.template getProp<1>(key)[i];

					//Check that the forces are (almost) equal
					float per = 0.1;
					if (a1 != 0.0)
						per = fabs(0.1*avg.get(i)/a1);

					BOOST_REQUIRE_CLOSE(a1,a2,per);
				}

				++it_v;
			}
		}
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_sym_cell_list_test )
{
	long int k = 4096*create_vcluster().getProcessingUnits();

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing vector symmetric cell-list k<=",k);

	// 3D test
	for ( ; k > 8*big_step ; k-= (k > 2*big_step)?big_step:small_step )
	{
		double r_cut = 0.1;

		// domain
		Box<3,double> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// ghost, big enough to contain the interaction radius
		Ghost<3,float> ghost(r_cut);

		vector_dist<3,double, aggregate<double> > vd(k,box,bc,ghost);

		{
		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			vd.getPos(p)[0] = (double)rand()/RAND_MAX;
			vd.getPos(p)[1] = (double)rand()/RAND_MAX;
			vd.getPos(p)[2] = (double)rand()/RAND_MAX;

			++it;
		}
		}

		vd.map();
		vd.ghost_get<>();

		// Get the Cell list structure
		auto NN = vd.getCellList(r_cut);

		// Get an iterator over particles
		auto it2 = vd.getDomainAndGhostIterator();

		openfpm::vector<openfpm::vector<size_t>> idx;
		idx.resize(vd.size_local_with_ghost());

		/////////// SYMMETRIC CASE CELL-LIST ////////

		// For each particle ...
		while (it2.isNext())
		{
			// ... p
			auto p = it2.get();

			// Get the position of the particle p
			Point<3,double> xp = vd.getPos(p);

			// Get an iterator over the neighborhood of the particle p symmetric
			auto NpSym = NN.template getNNIteratorSym<NO_CHECK>(NN.getCell(vd.getPos(p)),p.getKey());

			// For each neighborhood of the particle p
			while (NpSym.isNext())
			{
				// Neighborhood particle q
				auto q = NpSym.get();

				// if p == q skip this particle
				if (q == p.getKey() )	{++NpSym; continue;};

				// Get position of the particle q
				Point<3,double> xq = vd.getPos(q);

				// take the normalized direction
				double rn = norm2(xp - xq);

				// potential energy (using pow is slower)
				vd.getProp<0>(p) += rn;
				vd.getProp<0>(q) += rn;

				idx.get(p.getKey()).add(q);
				idx.get(q).add(p.getKey());


				// Next neighborhood
				++NpSym;
			}

			// Next Particle
			++it2;
		}

		/////////////// NON SYMMETRIC CASE ////////////////////////

		openfpm::vector<openfpm::vector<size_t>> idx2;
		idx2.resize(vd.size_local());

		auto it = vd.getDomainIterator();

		// For each particle ...
		while (it.isNext())
		{
			// ... p
			auto p = it.get();

			// Get the position of the particle p
			Point<3,double> xp = vd.getPos(p);

			// Get an iterator over the neighborhood of the particle p
			auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(p)));

			double Ep = 0.0;

			// For each neighborhood of the particle p
			while (Np.isNext())
			{
				// Neighborhood particle q
				auto q = Np.get();

				// if p == q skip this particle
				if (q == p.getKey())	{++Np; continue;};

				// Get position of the particle q
				Point<3,double> xq = vd.getPos(q);

				// take the normalized direction
				double rn = norm2(xp - xq);

				idx2.get(p.getKey()).add(q);

				// potential energy (using pow is slower)
				Ep += rn;

				// Next neighborhood
				++Np;
			}

			idx.get(p.getKey()).sort();
			idx2.get(p.getKey()).sort();

			bool ret = true;

			for (size_t i = 0 ; i < idx.get(p.getKey()).size() ; i++)
				ret &= idx.get(p.getKey()).get(i) == idx2.get(p.getKey()).get(i);

			BOOST_REQUIRE_EQUAL(ret,true);

			BOOST_REQUIRE_CLOSE(Ep,vd.getProp<0>(p),0.01);

			// Next Particle
			++it;
		}
	}
}


BOOST_AUTO_TEST_CASE( vector_dist_sym_verlet_list_test )
{
	if (create_vcluster().getProcessingUnits() > 12)
		return;

	long int k = 4096*create_vcluster().getProcessingUnits();

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	print_test( "Testing vector symmetric verlet-list k<=",k);

	// 3D test
	for ( ; k > 8*big_step ; k-= (k > 2*big_step)?big_step:small_step )
	{
		double r_cut = 0.1;

		// domain
		Box<3,double> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// ghost, big enough to contain the interaction radius
		Ghost<3,float> ghost(r_cut);

		vector_dist<3,double, aggregate<double> > vd(k,box,bc,ghost);

		{
		auto it = vd.getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();

			vd.getPos(p)[0] = (double)rand()/RAND_MAX;
			vd.getPos(p)[1] = (double)rand()/RAND_MAX;
			vd.getPos(p)[2] = (double)rand()/RAND_MAX;

			++it;
		}
		}

		vd.map();
		vd.ghost_get<>();

		// Get the Cell list structure
		auto NN = vd.getVerlet(r_cut,VL_SYMMETRIC);
		auto NNc = vd.getVerlet(r_cut);

		// Get an iterator over particles
		auto it2 = vd.getDomainAndGhostIterator();

		openfpm::vector<openfpm::vector<size_t>> idx;
		idx.resize(vd.size_local_with_ghost());

		/////////// SYMMETRIC CASE VERLET-LIST ////////

		// For each particle ...
		while (it2.isNext())
		{
			// ... p
			auto p = it2.get();

			// Get the position of the particle p
			Point<3,double> xp = vd.getPos(p);

			// Get an iterator over the neighborhood of the particle p symmetric
			auto NpSym = NN.template getNNIterator<NO_CHECK>(p.getKey());

			// For each neighborhood of the particle p
			while (NpSym.isNext())
			{
				// Neighborhood particle q
				auto q = NpSym.get();

				// if p == q skip this particle
				if (q == p.getKey())	{++NpSym; continue;};

				// Get position of the particle q
				Point<3,double> xq = vd.getPos(q);

				// take the normalized direction
				double rn = norm2(xp - xq);

				// potential energy (using pow is slower)
				vd.getProp<0>(p) += rn;
				vd.getProp<0>(q) += rn;

				idx.get(p.getKey()).add(q);
				idx.get(q).add(p.getKey());


				// Next neighborhood
				++NpSym;
			}

			// Next Particle
			++it2;
		}

		/////////////// NON SYMMETRIC CASE VERLET-LIST ////////////////////////

		openfpm::vector<openfpm::vector<size_t>> idx2;
		idx2.resize(vd.size_local());

		auto it = vd.getDomainIterator();

		// For each particle ...
		while (it.isNext())
		{
			// ... p
			auto p = it.get();

			// Get the position of the particle p
			Point<3,double> xp = vd.getPos(p);

			// Get an iterator over the neighborhood of the particle p
			auto Np = NNc.template getNNIterator<NO_CHECK>(p.getKey());

			double Ep = 0.0;

			// For each neighborhood of the particle p
			while (Np.isNext())
			{
				// Neighborhood particle q
				auto q = Np.get();

				// if p == q skip this particle
				if (q == p.getKey())	{++Np; continue;};

				// Get position of the particle q
				Point<3,double> xq = vd.getPos(q);

				// take the normalized direction
				double rn = norm2(xp - xq);

				idx2.get(p.getKey()).add(q);

				// potential energy (using pow is slower)
				Ep += rn;

				// Next neighborhood
				++Np;
			}

			idx.get(p.getKey()).sort();
			idx2.get(p.getKey()).sort();

			bool ret = true;

			BOOST_REQUIRE_EQUAL(idx.get(p.getKey()).size(),idx2.get(p.getKey()).size());

			for (size_t i = 0 ; i < idx.get(p.getKey()).size() ; i++)
			{
				ret &= idx.get(p.getKey()).get(i) == idx2.get(p.getKey()).get(i);
			}

			BOOST_REQUIRE_EQUAL(ret,true);

			BOOST_REQUIRE_CLOSE(Ep,vd.getProp<0>(p),0.01);

			// Next Particle
			++it;
		}
	}
}


#endif /* SRC_VECTOR_VECTOR_DIST_CELL_LIST_TESTS_HPP_ */
