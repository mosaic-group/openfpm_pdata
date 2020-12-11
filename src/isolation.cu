#include <iostream>
#include <thread>

size_t debug_tot_call = 0;

#define PRINT_STACKTRACE
#define CHECKFOR_POSNAN
#define CHECKFOR_POSINF
#define CHECKFOR_PROPNAN
#define CHECKFOR_PROPINF

#define NO_WARNING
#include "Graph/CartesianGraphFactory.hpp"

void timeout_cycle()
{
	// 6 seconds
	std::this_thread::sleep_for (std::chrono::seconds(900));

	std::cout << "Time Out" << std::endl;
	std::exit(1);
}


#define BOOST_DISABLE_ASSERTS


#include "config.h"
#undef VERSION

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include <Vector/vector_dist.hpp>
#include "Vector/tests/vector_dist_util_unit_tests.hpp"

// initialization function:
bool init_unit_test()
{
//  std::thread to (timeout_cycle);
//  to.detach();
  return true;
}

// entry point
int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

  return boost::unit_test::unit_test_main( &init_unit_test, argc, argv );
}



BOOST_AUTO_TEST_CASE( vector_dist_ghost_put_gpu )
{
	

	Vcluster<> & v_cl = create_vcluster();

	long int k = 25*25*25*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	if (v_cl.getProcessingUnits() > 48)
		return;

	BOOST_TEST_CHECKPOINT( "Testing 3D periodic ghost put k=" << k );

	long int big_step = k / 30;
	big_step = (big_step == 0)?1:big_step;
	long int small_step = 21;

	// 3D test
	for ( ; k >= 2 ; k-= (k > 2*big_step)?big_step:small_step )
	{
		float r_cut = 1.3 / k;
		float r_g = 1.5 / k;

		Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

		// Boundary conditions
		size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

		// ghost
		Ghost<3,float> ghost(r_g);

		typedef  aggregate<float,float,float> part_prop;

		// Distributed vector
		vector_dist_gpu<3,float, part_prop > vd(0,box,bc,ghost);

		auto it = vd.getGridIterator({(size_t)k,(size_t)k,(size_t)k});

		while (it.isNext())
		{
			auto key = it.get();

			vd.add();

			vd.getLastPosWrite()[0] = key.get(0)*it.getSpacing(0);
			vd.getLastPosWrite()[1] = key.get(1)*it.getSpacing(1);
			vd.getLastPosWrite()[2] = key.get(2)*it.getSpacing(2);

			// Fill some properties randomly

			vd.getLastPropWrite<0>() = 0.0;

			vd.getLastPropWrite<2>() = 0.0;

			++it;
		}

		vd.map();

		vd.hostToDevicePos();
		vd.template hostToDeviceProp<0,2>();
		// sync the ghost
		vd.ghost_get<0,2>(RUN_ON_DEVICE);
		vd.template deviceToHostProp<0,2>();
		vd.deviceToHostPos();

		{
			auto NN = vd.getCellList(r_cut);
			float a = 1.0f*k*k;

			// run trough all the particles + ghost

			auto it2 = vd.getDomainIterator();

			while (it2.isNext())
			{
				// particle p
				auto p = it2.get();
				Point<3,float> xp = vd.getPos(p);

				// Get an iterator over the neighborhood particles of p
				auto Np = NN.getNNIterator<NO_CHECK>(NN.getCell(xp));

				// For each neighborhood particle ...
				while (Np.isNext())
				{
					auto q = Np.get();
					Point<3,float> xq = vd.getPosRead(q);

					float dist = xp.distance(xq);

					if (dist < r_cut)
					{
						vd.getPropWrite<0>(q) += a*(-dist*dist+r_cut*r_cut);
						vd.getPropWrite<2>(q) += a*(-dist*dist+r_cut*r_cut);
					}

					++Np;
				}

				++it2;
			}

			vd.hostToDevicePos();
			vd.template hostToDeviceProp<0,2>();
			vd.template ghost_put<add_atomic_,0,2>(RUN_ON_DEVICE);
			vd.template deviceToHostProp<0,2>();
			vd.deviceToHostPos();

			bool ret = true;
			auto it3 = vd.getDomainIterator();

			float constant = vd.getProp<0>(it3.get());
			float eps = 0.001;

			while (it3.isNext())
			{
				float constant2 = vd.getProp<0>(it3.get());
				float constant3 = vd.getProp<2>(it3.get());
				if (fabs(constant - constant2)/constant > eps || fabs(constant - constant3)/constant > eps)
				{
					Point<3,float> p = vd.getPosRead(it3.get());

					std::cout << p.toString() << "    " <<  constant2 << "/" << constant << "/" << constant3 << "    " << v_cl.getProcessUnitID() << std::endl;
					ret = false;
					break;
				}

				++it3;
			}
			BOOST_REQUIRE_EQUAL(ret,true);
		}

		auto itp = vd.getDomainAndGhostIterator();
		while (itp.isNext())
		{
			auto key = itp.get();

			vd.getPropWrite<0>(key) = 0.0;
			vd.getPropWrite<2>(key) = 0.0;

			++itp;
		}

		{
			auto NN = vd.getCellList(r_cut);
			float a = 1.0f*k*k;

			// run trough all the particles + ghost

			auto it2 = vd.getDomainIterator();

			while (it2.isNext())
			{
				// particle p
				auto p = it2.get();
				Point<3,float> xp = vd.getPosRead(p);

				// Get an iterator over the neighborhood particles of p
				auto Np = NN.getNNIterator<NO_CHECK>(NN.getCell(xp));

				// For each neighborhood particle ...
				while (Np.isNext())
				{
					auto q = Np.get();
					Point<3,float> xq = vd.getPosRead(q);

					float dist = xp.distance(xq);

					if (dist < r_cut)
					{
						vd.getPropWrite<0>(q) += a*(-dist*dist+r_cut*r_cut);
						vd.getPropWrite<2>(q) += a*(-dist*dist+r_cut*r_cut);
					}

					++Np;
				}

				++it2;
			}

			vd.hostToDevicePos();
			vd.template hostToDeviceProp<0,2>();
			vd.template ghost_put<add_atomic_,0>(RUN_ON_DEVICE);
			vd.template ghost_put<add_atomic_,2>(RUN_ON_DEVICE);
			vd.template deviceToHostProp<0,2>();
			vd.deviceToHostPos();

			bool ret = true;
			auto it3 = vd.getDomainIterator();

			float constant = vd.getPropRead<0>(it3.get());
			float eps = 0.001;

			while (it3.isNext())
			{
				float constant2 = vd.getPropRead<0>(it3.get());
				float constant3 = vd.getPropRead<0>(it3.get());
				if (fabs(constant - constant2)/constant > eps || fabs(constant - constant3)/constant > eps)
				{
					Point<3,float> p = vd.getPosRead(it3.get());

					std::cout << p.toString() << "    " <<  constant2 << "/" << constant << "/" << constant3 << "  " << it3.get().getKey() << "    " << v_cl.getProcessUnitID() << std::endl;
					ret = false;
					break;
				}

				++it3;
			}
			BOOST_REQUIRE_EQUAL(ret,true);
		}
	}

	openfpm_finalize();
}