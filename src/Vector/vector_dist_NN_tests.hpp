/*
 * vector_dist_NN_tests.hpp
 *
 *  Created on: Aug 20, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_NN_TESTS_HPP_
#define SRC_VECTOR_VECTOR_DIST_NN_TESTS_HPP_


BOOST_AUTO_TEST_CASE( vector_dist_full_NN )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 12)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 750 * v_cl.getProcessingUnits();

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test("Testing 3D full NN search k=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D full NN search k=" << k );

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	for (double r_cut = 0.1 ; r_cut < 1.0; r_cut += 0.1)
	{
		// ghost
		Ghost<3,float> ghost(r_cut*1.001);

		typedef  aggregate<float> part_prop;

		// Distributed vector
		vector_dist<3,float, part_prop > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			// Fill some properties randomly

			vd.getProp<0>(key) = 0.0;

			++it;
		}

		vd.map();

		// sync the ghost
		vd.ghost_get<0>();

		openfpm::vector<openfpm::vector<size_t>> list_idx;
		openfpm::vector<openfpm::vector<size_t>> list_idx2;

		list_idx.resize(vd.size_local());
		list_idx2.resize(vd.size_local());

		for (size_t i = 0 ; i < vd.size_local() ; i++)
		{
			Point<3,float> p = vd.getPos(i);

			for (size_t j = 0 ; j < vd.size_local_with_ghost(); j++)
			{
				Point<3,float> q = vd.getPos(j);

				if (p.distance2(q) < r_cut * r_cut)
					list_idx.get(i).add(j);
			}

			list_idx.get(i).sort();
		}

		auto NN = vd.getCellList(r_cut);

		it = vd.getDomainIterator();

		while (it.isNext())
		{
			Point<3,float> xp = vd.getPos(it.get());
			auto Np = NN.getNNIterator<NO_CHECK>(NN.getCell(xp));

			while (Np.isNext())
			{
				auto q = Np.get();

				Point<3,float> xq = vd.getPos(q);

				if (xp.distance2(xq) < r_cut * r_cut)
					list_idx2.get(it.get().getKey()).add(q);

				++Np;
			}

			list_idx2.get(it.get().getKey()).sort();

			++it;
		}

		BOOST_REQUIRE_EQUAL(list_idx.size(),list_idx2.size());

		bool ret;
		for (size_t i = 0 ; i < list_idx.size() ; i++)
		{
			BOOST_REQUIRE_EQUAL(list_idx.get(i).size(),list_idx2.get(i).size());

			ret = true;
			for (size_t j = 0 ; j < list_idx.get(i).size() ; j++)
				ret &= list_idx.get(i).get(j) == list_idx2.get(i).get(j);

			BOOST_REQUIRE_EQUAL(ret,true);
		}

		///////////////////////////////////

		auto NNv = vd.getVerlet(r_cut*1.0001);

		it = vd.getDomainIterator();

		while (it.isNext())
		{
			Point<3,float> xp = vd.getPos(it.get());
			auto Np = NNv.getNNIterator<NO_CHECK>(it.get().getKey());

			list_idx2.get(it.get().getKey()).clear();

			while (Np.isNext())
			{
				auto q = Np.get();

				Point<3,float> xq = vd.getPos(q);

				if (xp.distance2(xq) < r_cut * r_cut)
					list_idx2.get(it.get().getKey()).add(q);

				++Np;
			}

			list_idx2.get(it.get().getKey()).sort();

			++it;
		}

		BOOST_REQUIRE_EQUAL(list_idx.size(),list_idx2.size());

		for (size_t i = 0 ; i < list_idx.size() ; i++)
		{
			BOOST_REQUIRE_EQUAL(list_idx.get(i).size(),list_idx2.get(i).size());

			ret = true;
			for (size_t j = 0 ; j < list_idx.get(i).size() ; j++)
				ret &= list_idx.get(i).get(j) == list_idx2.get(i).get(j);

			BOOST_REQUIRE_EQUAL(ret,true);
		}

		//////////////////////////////////////////

		NNv.clear();
		vd.updateVerlet(NNv,r_cut*1.0001);

		it = vd.getDomainIterator();

		while (it.isNext())
		{
			Point<3,float> xp = vd.getPos(it.get());
			auto Np = NNv.getNNIterator<NO_CHECK>(it.get().getKey());

			list_idx2.get(it.get().getKey()).clear();

			while (Np.isNext())
			{
				auto q = Np.get();

				Point<3,float> xq = vd.getPos(q);

				if (xp.distance2(xq) < r_cut * r_cut)
					list_idx2.get(it.get().getKey()).add(q);

				++Np;
			}

			list_idx2.get(it.get().getKey()).sort();

			++it;
		}

		BOOST_REQUIRE_EQUAL(list_idx.size(),list_idx2.size());

		for (size_t i = 0 ; i < list_idx.size() ; i++)
		{
			BOOST_REQUIRE_EQUAL(list_idx.get(i).size(),list_idx2.get(i).size());

			ret = true;
			for (size_t j = 0 ; j < list_idx.get(i).size() ; j++)
				ret &= list_idx.get(i).get(j) == list_idx2.get(i).get(j);

			BOOST_REQUIRE_EQUAL(ret,true);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_particle_iteration )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 12)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 750 * v_cl.getProcessingUnits();

	print_test("Testing 3D particle cell iterator=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D full NN search k=" << k );

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 0.01;

	// ghost
	Ghost<3,float> ghost(r_cut);

	typedef  aggregate<float> part_prop;

	// Distributed vector
	vector_dist<3,float, part_prop > vd(k,box,bc,ghost,BIND_DEC_TO_GHOST);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			// Fill some properties randomly

			vd.getProp<0>(key) = 0.0;

			++it;
		}

	vd.map();

	// sync the ghost
	vd.ghost_get<0>();

	openfpm::vector<size_t> ids;
	ids.resize(vd.size_local());

	auto NN = vd.getCellListSym(r_cut);

	auto it_pcell = vd.getDomainIteratorCells(NN);

	size_t count = 0;
	while (it_pcell.isNext())
	{
		count++;

		size_t id = it_pcell.get();
		ids.get(id) = 1;

		BOOST_REQUIRE(id < vd.size_local());

		++it_pcell;
	}

	v_cl.sum(count);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL((long int)count,k);
}

#endif /* SRC_VECTOR_VECTOR_DIST_NN_TESTS_HPP_ */
