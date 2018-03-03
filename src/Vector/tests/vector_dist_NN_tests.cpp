/*
 * vector_dist_NN_tests.hpp
 *
 *  Created on: Aug 20, 2016
 *      Author: i-bird
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "VCluster/VCluster.hpp"
#include "Vector/vector_dist.hpp"

extern void print_test_v(std::string test, size_t sz);

template<typename VerletList>
void test_full_nn(long int k)
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 12)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v("Testing 3D full NN search k=",k);
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

		auto NNv = vd.template getVerlet<VerletList>(r_cut*1.0001);

		it = vd.getDomainIterator();

		while (it.isNext())
		{
			Point<3,float> xp = vd.getPos(it.get());
			auto Np = NNv.template getNNIterator<NO_CHECK>(it.get().getKey());

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
			auto Np = NNv.template getNNIterator<NO_CHECK>(it.get().getKey());

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

BOOST_AUTO_TEST_CASE( vector_dist_full_NN )
{
	auto & v_cl = create_vcluster();

#ifdef TEST_COVERAGE_MODE
    long int k = 50 * v_cl.getProcessingUnits();
#else
    long int k = 750 * v_cl.getProcessingUnits();
#endif

	test_full_nn<VERLET_MEMFAST(3,float)>(k);

	k /= 2;
	test_full_nn<VERLET_MEMBAL(3,float)>(k);
	k /= 2;
	test_full_nn<VERLET_MEMMW(3,float)>(k);
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

	print_test_v("Testing 3D particle cell iterator=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D full NN search k=" << k );

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 0.1;

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


BOOST_AUTO_TEST_CASE( vector_dist_particle_NN_MP_iteration )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 24)
	{return;}

	float L = 1000.0;

    // set the seed
	// create the random generator engine
	std::srand(0);
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(-L,L);

    long int k = 4096 * v_cl.getProcessingUnits();

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v("Testing 3D periodic vector symmetric cell-list k=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D periodic vector symmetric cell-list k=" << k );

	Box<3,float> box({-L,-L,-L},{L,L,L});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 100.0;

	// ghost
	Ghost<3,float> ghost(r_cut);

	// Point and global id
	struct point_and_gid
	{
		size_t id;
		Point<3,float> xq;

		bool operator<(const struct point_and_gid & pag) const
		{
			return (id < pag.id);
		}
	};

	typedef  aggregate<size_t,size_t,size_t,openfpm::vector<point_and_gid>,openfpm::vector<point_and_gid>> part_prop;

	// Distributed vector
	vector_dist<3,float, part_prop > vd(k,box,bc,ghost,BIND_DEC_TO_GHOST);
	size_t start = vd.init_size_accum(k);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPosWrite(key)[0] = ud(eg);
		vd.getPosWrite(key)[1] = ud(eg);
		vd.getPosWrite(key)[2] = ud(eg);

		// Fill some properties randomly

		vd.getPropWrite<0>(key) = 0;
		vd.getPropWrite<1>(key) = 0;
		vd.getPropWrite<2>(key) = key.getKey() + start;

		++it;
	}

	vd.map();

	// sync the ghost
	vd.ghost_get<0,2>();

	auto NN = vd.getCellList(r_cut);
	auto p_it = vd.getDomainIterator();

	while (p_it.isNext())
	{
		auto p = p_it.get();

		Point<3,float> xp = vd.getPosRead(p);

		auto Np = NN.getNNIterator(NN.getCell(xp));

		while (Np.isNext())
		{
			auto q = Np.get();

			if (p.getKey() == q)
			{
				++Np;
				continue;
			}

			// repulsive

			Point<3,float> xq = vd.getPosRead(q);
			Point<3,float> f = (xp - xq);

			float distance = f.norm();

			// Particle should be inside 2 * r_cut range

			if (distance < r_cut )
			{
				vd.getPropWrite<0>(p)++;
				vd.getPropWrite<3>(p).add();
				vd.getPropWrite<3>(p).last().xq = xq;
				vd.getPropWrite<3>(p).last().id = vd.getPropWrite<2>(q);
			}

			++Np;
		}

		++p_it;
	}

	// We now divide the particles on 4 phases

	openfpm::vector<vector_dist<3,float, part_prop >> phases;
	phases.add( vector_dist<3,float, part_prop >(vd.getDecomposition(),0));
	phases.add( vector_dist<3,float, part_prop >(phases.get(0).getDecomposition(),0));
	phases.add( vector_dist<3,float, part_prop >(phases.get(0).getDecomposition(),0));
	phases.add( vector_dist<3,float, part_prop >(phases.get(0).getDecomposition(),0));

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		if (p.getKey() % 4 == 0)
		{
			phases.get(0).add();
			phases.get(0).getLastPos()[0] = vd.getPos(p)[0];
			phases.get(0).getLastPos()[1] = vd.getPos(p)[1];
			phases.get(0).getLastPos()[2] = vd.getPos(p)[2];

			phases.get(0).template getLastProp<2>() = vd.template getProp<2>(p);
		}
		else if (p.getKey() % 4 == 1)
		{
			phases.get(1).add();
			phases.get(1).getLastPos()[0] = vd.getPos(p)[0];
			phases.get(1).getLastPos()[1] = vd.getPos(p)[1];
			phases.get(1).getLastPos()[2] = vd.getPos(p)[2];

			phases.get(1).template getLastProp<2>() = vd.template getProp<2>(p);
		}
		else if (p.getKey() % 4 == 2)
		{
			phases.get(2).add();
			phases.get(2).getLastPos()[0] = vd.getPos(p)[0];
			phases.get(2).getLastPos()[1] = vd.getPos(p)[1];
			phases.get(2).getLastPos()[2] = vd.getPos(p)[2];

			phases.get(2).template getLastProp<2>() = vd.template getProp<2>(p);
		}
		else
		{
			phases.get(3).add();
			phases.get(3).getLastPos()[0] = vd.getPos(p)[0];
			phases.get(3).getLastPos()[1] = vd.getPos(p)[1];
			phases.get(3).getLastPos()[2] = vd.getPos(p)[2];

			phases.get(3).template getLastProp<2>() = vd.template getProp<2>(p);
		}

		++it2;
	}

	// now we get all the Cell-lists

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		phases.get(i).ghost_get<0,2>();
	}

	openfpm::vector<CellList<3, float, Mem_fast<>, shift<3, float> >> NN_ptr;

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		NN_ptr.add(phases.get(i).getCellListSym(r_cut));
	}

	// We now interact all the phases

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		for (size_t j = 0 ; j < phases.size() ; j++)
		{
			auto p_it2 = phases.get(i).getDomainIterator();

			while (p_it2.isNext())
			{
				auto p = p_it2.get();

				Point<3,float> xp = phases.get(i).getPosRead(p);

				auto Np = NN_ptr.get(j).getNNIteratorSymMP<NO_CHECK>(NN_ptr.get(j).getCell(xp),p.getKey(),phases.get(i).getPosVector(),phases.get(j).getPosVector());

				while (Np.isNext())
				{
					auto q = Np.get();

					if (p.getKey() == q && i == j)
					{
						++Np;
						continue;
					}

					// repulsive

					Point<3,float> xq = phases.get(j).getPosRead(q);
					Point<3,float> f = (xp - xq);

					float distance = f.norm();

					// Particle should be inside r_cut range

					if (distance < r_cut )
					{
						phases.get(i).getPropWrite<1>(p)++;
						phases.get(j).getPropWrite<1>(q)++;

						phases.get(i).getPropWrite<4>(p).add();
						phases.get(j).getPropWrite<4>(q).add();

						phases.get(i).getPropWrite<4>(p).last().xq = xq;
						phases.get(j).getPropWrite<4>(q).last().xq = xp;
						phases.get(i).getPropWrite<4>(p).last().id = phases.get(j).getProp<2>(q);
						phases.get(j).getPropWrite<4>(q).last().id = phases.get(i).getProp<2>(p);
					}

					++Np;
				}

				++p_it2;
			}
		}
	}

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		phases.get(i).ghost_put<add_,1>();
		phases.get(i).ghost_put<merge_,4>();
	}

	auto p_it3 = vd.getDomainIterator();

	bool ret = true;
	while (p_it3.isNext())
	{
		auto p = p_it3.get();

		int ph;

		if (p.getKey() % 4 == 0)
		{ph = 0;}
		else if (p.getKey() % 4 == 1)
		{ph = 1;}
		else if (p.getKey() % 4 == 2)
		{ph = 2;}
		else
		{ph = 3;}

		size_t pah = p.getKey()/4;
		ret &= phases.get(ph).getPropRead<1>(pah) == vd.getPropRead<0>(p);

		vd.getPropRead<3>(p).sort();
		phases.get(ph).getPropRead<4>(pah).sort();

		ret &= vd.getPropRead<3>(p).size() == phases.get(ph).getPropRead<4>(pah).size();

		for (size_t i = 0 ; i < vd.getPropRead<3>(p).size() ; i++)
			ret &= vd.getPropRead<3>(p).get(i).id == phases.get(ph).getPropRead<4>(pah).get(i).id;

		if (ret == false)
		{
			std::cout << "Error on particle: " << vd.getPropRead<2>(p) << "   " << v_cl.rank() << std::endl;

			std::cout << vd.getPropRead<3>(p).size() << "   " << phases.get(ph).getPropRead<4>(pah).size() << "  " << v_cl.rank() << std::endl;

			for (size_t i = 0 ; i < vd.getPropRead<3>(p).size() ; i++)
				std::cout << vd.getPropRead<3>(p).get(i).id << "    " << phases.get(ph).getPropRead<4>(pah).get(i).id << "  " << v_cl.rank() << std::endl;

			std::cout << phases.get(ph).getPropRead<1>(pah) << "  A  " << vd.getPropRead<0>(p) << std::endl;

			break;
		}

		++p_it3;
	}

	BOOST_REQUIRE_EQUAL(ret,true);
}
