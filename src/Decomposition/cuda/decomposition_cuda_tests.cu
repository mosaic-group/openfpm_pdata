#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "VCluster/VCluster.hpp"
#include "Decomposition/CartDecomposition.hpp"

#define SUB_UNIT_FACTOR 1024

template<typename dec_type>
__global__ void test_proc_idbc(Point<3,double> p1 ,Point<3,double> p2 , dec_type dec, unsigned int * pr_id)
{
	pr_id[0] = dec.processorIDBC(p1);
	pr_id[1] = dec.processorIDBC(p2);
}

template<typename dec_type>
__global__ void test_ghost_n(Point<3,double> p1 ,Point<3,double> p2 , dec_type dec, unsigned int * ng_id)
{
	ng_id[0] = dec.ghost_processorID_N(p1);
	ng_id[1] = dec.ghost_processorID_N(p2);
}

template<typename dec_type, typename output_type>
__global__ void test_ghost(Point<3,double> p1 ,Point<3,double> p2 , dec_type dec, unsigned int * ng_id , output_type g_id)
{
	for (unsigned int i = 0 ; i < ng_id[0] ; i++)
	{
		dec.ghost_processor_ID(p1,g_id,0,i);
	}

	for (unsigned int i = 0 ; i < ng_id[1] ; i++)
	{
		dec.ghost_processor_ID(p2,g_id,ng_id[0],i);
	}
}

BOOST_AUTO_TEST_SUITE( decomposition_to_gpu_test )

BOOST_AUTO_TEST_CASE( CartDecomposition_check_cross_consistency_between_proc_idbc_and_ghost2_gpu )
{
	// Vcluster
	Vcluster<> & vcl = create_vcluster();

	CartDecomposition<3, double, CudaMemory,memory_traits_inte> dec(vcl);

	size_t bc[3] = {PERIODIC,PERIODIC,PERIODIC};

	// Physical domain
	Box<3, double> box( { -0.01, -0.01, 0.0 }, { 0.01, 0.01, 0.003 });

	Ghost<3,double> g(0.0015);

	dec.setGoodParameters(box, bc, g, 512);

	dec.decompose();

	// Now we check the point

	for (size_t j = 0 ; j < 3 ; j++ )
	{
		for (size_t i = 0 ; i < dec.getNSubDomain() ; i++)
		{
			Point<3,double> p1;
			Point<3,double> p2;

			p1.get(0) = SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(0);
			p1.get(1) = SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(1);
			p1.get(2) = SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(2);

			p2 = p1;

//			p2.get(j) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(j),-1.0);

			auto gpudec = dec.toKernel();

			CudaMemory mem;
			mem.allocate(2*sizeof(unsigned int));

			test_proc_idbc<decltype(gpudec)><<<1,1>>>(p1,p2,gpudec,(unsigned int *)mem.getDevicePointer());

			mem.deviceToHost();

			BOOST_REQUIRE(((unsigned int *)mem.getPointer())[0] < vcl.size());
			BOOST_REQUIRE(((unsigned int *)mem.getPointer())[1] < vcl.size());

			CudaMemory mem2;
			mem2.allocate(2*sizeof(unsigned int));
			test_ghost_n<decltype(gpudec)><<<1,1>>>(p1,p2,gpudec,(unsigned int *)mem2.getDevicePointer());

			mem2.deviceToHost();

			unsigned int tot = ((unsigned int *)mem2.getPointer())[0] + ((unsigned int *)mem2.getPointer())[1];

			openfpm::vector_gpu<aggregate<int,int>> vd;
			vd.resize(tot);
			test_ghost<decltype(gpudec),decltype(vd.toKernel())><<<1,1>>>(p1,p2,gpudec,(unsigned int *)mem2.getDevicePointer(),vd.toKernel());

			if (((unsigned int *)mem.getPointer())[0] != ((unsigned int *)mem.getPointer())[1])
			{
				if (vcl.rank() == ((unsigned int *)mem.getPointer())[1] )
				{
					BOOST_REQUIRE(((unsigned int *)mem2.getPointer())[1] != 0);
					BOOST_REQUIRE(((unsigned int *)mem2.getPointer())[0] == 0);
				}

				if (vcl.rank() == ((unsigned int *)mem.getPointer())[0])
				{
					BOOST_REQUIRE(((unsigned int *)mem2.getPointer())[1] == 0 );
					BOOST_REQUIRE(((unsigned int *)mem2.getPointer())[0] != 0 );
				}
			}


			p1.get(0) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getHigh(0),SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(0));
			p1.get(1) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getHigh(1),SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(1));
			p1.get(2) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getHigh(2),SpaceBox<3,double>(dec.getSubDomains().get(i)).getLow(2));

			p2 = p1;

			p2.get(j) = std::nextafter(SpaceBox<3,double>(dec.getSubDomains().get(i)).getHigh(j),1.0);

			test_proc_idbc<decltype(gpudec)><<<1,1>>>(p1,p2,gpudec,(unsigned int *)mem.getDevicePointer());

			mem.deviceToHost();

			BOOST_REQUIRE(((unsigned int *)mem.getPointer())[0] < vcl.size());
			BOOST_REQUIRE(((unsigned int *)mem.getPointer())[1] < vcl.size());

			mem2.allocate(2*sizeof(unsigned int));
			test_ghost_n<decltype(gpudec)><<<1,1>>>(p1,p2,gpudec,(unsigned int *)mem2.getDevicePointer());

			mem2.deviceToHost();

			tot = ((unsigned int *)mem2.getPointer())[0] + ((unsigned int *)mem2.getPointer())[1];

			vd.resize(tot);
			test_ghost<decltype(gpudec),decltype(vd.toKernel())><<<1,1>>>(p1,p2,gpudec,(unsigned int *)mem2.getDevicePointer(),vd.toKernel());

			if (((unsigned int *)mem.getPointer())[0] != ((unsigned int *)mem.getPointer())[1])
			{
				if (vcl.rank() == ((unsigned int *)mem.getPointer())[1])
				{
					BOOST_REQUIRE(((unsigned int *)mem2.getPointer())[1] != 0);
					BOOST_REQUIRE(((unsigned int *)mem2.getPointer())[0] == 0);
				}

				if (vcl.rank() == ((unsigned int *)mem.getPointer())[0])
				{
					BOOST_REQUIRE(((unsigned int *)mem2.getPointer())[1] == 0 );
					BOOST_REQUIRE(((unsigned int *)mem2.getPointer())[0] != 0 );
				}
			}

		}
	}
}

BOOST_AUTO_TEST_SUITE_END()
