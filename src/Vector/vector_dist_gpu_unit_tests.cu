
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "VCluster/VCluster.hpp"
#include <Vector/vector_dist.hpp>


BOOST_AUTO_TEST_SUITE( vector_dist_gpu_test )

void print_test(std::string test, size_t sz)
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << test << " " << sz << "\n";
}


__global__  void initialize_props(vector_dist_ker<3, float, aggregate<float, float [3], float[3]>> vd)
{
	auto p = GET_PARTICLE(vd);

	vd.template getProp<0>(p) = vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2];

	vd.template getProp<1>(p)[0] = vd.getPos(p)[0] + vd.getPos(p)[1];
	vd.template getProp<1>(p)[1] = vd.getPos(p)[0] + vd.getPos(p)[2];
	vd.template getProp<1>(p)[2] = vd.getPos(p)[1] + vd.getPos(p)[2];
}

template<typename CellList_type>
__global__  void calculate_force(vector_dist_ker<3, float, aggregate<float, float[3], float [3]>> vd,
		                         vector_dist_ker<3, float, aggregate<float, float[3], float [3]>> vd_sort,
		                         CellList_type cl)
{
	auto p = GET_PARTICLE(vd);

	Point<3,float> xp = vd.getPos(p);

    auto it = cl.getNNIterator(cl.getCell(xp));

    auto cell = cl.getCell(xp);

    int s1 = cell.get(0);
    int s2 = cell.get(1);
    int s3 = cell.get(2);

    Point<3,float> force1({0.0,0.0,0.0});
    Point<3,float> force2({0.0,0.0,0.0});

    while (it.isNext())
    {
    	auto q1 = it.get();
    	auto q2 = it.get_orig();

    	if (q2 == p) {++it; continue;}

    	Point<3,float> xq_1 = vd_sort.getPos(q1);
    	Point<3,float> xq_2 = vd.getPos(q2);

    	Point<3,float> r1 = xq_1 - xp;
    	Point<3,float> r2 = xq_2 - xp;

    	// Normalize

    	r1 /= r1.norm();
    	r2 /= r2.norm();

    	force1 += vd_sort.template getProp<0>(q1)*r1;
    	force2 += vd.template getProp<0>(q2)*r2;

    	++it;
    }

    vd.template getProp<1>(p)[0] = force1.get(0);
    vd.template getProp<1>(p)[1] = force1.get(1);
    vd.template getProp<1>(p)[2] = force1.get(2);

    vd.template getProp<2>(p)[0] = force2.get(0);
    vd.template getProp<2>(p)[1] = force2.get(1);
    vd.template getProp<2>(p)[2] = force2.get(2);
}

BOOST_AUTO_TEST_CASE( vector_dist_gpu_test)
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 16)
	{return;}

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// set the ghost based on the radius cut off (make just a little bit smaller than the spacing)
	Ghost<3,float> g(0.1);

	// Boundary conditions
	size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	vector_dist_gpu<3,float,aggregate<float,float[3],float[3]>> vd(1000,domain,bc,g);

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0] = (float)rand() / RAND_MAX;
		vd.getPos(p)[1] = (float)rand() / RAND_MAX;
		vd.getPos(p)[2] = (float)rand() / RAND_MAX;

		++it;
	}

	// Ok we redistribute the particles
	vd.map();

	size_t size_l = vd.size_local();

	v_cl.sum(size_l);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(size_l,1000);


	auto & ct = vd.getDecomposition();

	bool noOut = true;
	size_t cnt = 0;

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		noOut &= ct.isLocal(vd.getPos(p));

		cnt++;
		++it2;
	}

	BOOST_REQUIRE_EQUAL(noOut,true);
	BOOST_REQUIRE_EQUAL(cnt,vd.size_local());

	vd.write("test_out_gpu");

	// now we offload all the properties

	auto it3 = vd.getDomainIteratorGPU();

	vd.hostToDevicePos();

	initialize_props<<<it3.wthr,it3.thr>>>(vd.toKernel());

	// now we check what we initialized

	vd.deviceToHostProp<0,1>();

	auto it4 = vd.getDomainIterator();

	while (it4.isNext())
	{
		auto p = it4.get();

		BOOST_REQUIRE_CLOSE(vd.template getProp<0>(p),vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2],0.01);

		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[0],vd.getPos(p)[0] + vd.getPos(p)[1],0.01);
		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[1],vd.getPos(p)[0] + vd.getPos(p)[2],0.01);
		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[2],vd.getPos(p)[1] + vd.getPos(p)[2],0.01);

		//std::cout << "PROP 0 " << vd.template getProp<0>(p) << "   " << vd.getPos(p)[0] + vd.getPos(p)[1] + vd.getPos(p)[2] << std::endl;

		++it4;
	}

	auto NN = vd.getCellListGPU(0.1);
	auto NN_cpu = vd.getCellList(0.1);

	auto it5 = vd.getDomainIteratorGPU();

	calculate_force<decltype(NN.toKernel())><<<it5.wthr,it5.thr>>>(vd.toKernel(),vd.toKernel_sorted(),NN.toKernel());

	vd.template deviceToHostProp<1,2>();

	auto it6 = vd.getDomainIterator();

	bool match = true;

	while (it6.isNext())
	{
		auto p = it6.get();

		Point<3,float> xp = vd.getPos(p);

		// Calculate on CPU

		Point<3,float> force({0.0,0.0,0.0});

		auto NNc = NN_cpu.getNNIterator(NN_cpu.getCell(xp));

		while (NNc.isNext())
		{
			auto q = NNc.get();

	    	if (q == p.getKey()) {++NNc; continue;}

	    	Point<3,float> xq_2 = vd.getPos(q);
	    	Point<3,float> r2 = xq_2 - xp;

	    	// Normalize

	    	r2 /= r2.norm();
	    	force += vd.template getProp<0>(q)*r2;

			++NNc;
		}

		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[0],vd.template getProp<2>(p)[0],0.001);
		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[1],vd.template getProp<2>(p)[1],0.001);
		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[2],vd.template getProp<2>(p)[2],0.001);

		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[0],force.get(0),0.01);
		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[1],force.get(1),0.01);
		BOOST_REQUIRE_CLOSE(vd.template getProp<1>(p)[2],force.get(2),0.01);

		++it6;
	}

}

BOOST_AUTO_TEST_SUITE_END()
