#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Point_test.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"

extern void print_test_v(std::string test, size_t sz);

BOOST_AUTO_TEST_SUITE( grid_dist_id_test )


BOOST_AUTO_TEST_CASE( grid_dist_id_gpu_test )
{
	// Test grid periodic

/*	Box<3,float> domain({-1.0,-1.0,-1.0},{1.0,1.0,1.0});

	Vcluster<> & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
	{return;}

	// grid size
	size_t sz[3];
	sz[0] = 32;
	sz[1] = 32;
	sz[2] = 32;

	// Ghost
	Ghost<3,long int> g(1);

	// periodicity
	periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

	// Distributed grid with id decomposition
    grid_dist_id_gpu<3, float, aggregate<float, float>> g_dist(sz,domain,g,pr);
   
	void * ptr = (double*)g_dist.get_loc_grid(0).toKernel().getPointer<0>();*/

/*	Box<3,size_t> box({1,1,1},{30,30,30});
    auto it = g_dist.getGridIterator(box.getKP1(),box.getKP2());

    float c = 5.0;

    typedef typename GetSetBlockType<decltype(g_dist)>::type BlockT;

	g_dist.setPoints(box.getKP1(),box.getKP2(),
			        [c] __device__ (BlockT & data, int i, int j, int k)
			        {
			        	data.template get<0>() = c + i*i + j*j + k*k;
			        }
                    );
  */                  
    
}


BOOST_AUTO_TEST_SUITE_END()
