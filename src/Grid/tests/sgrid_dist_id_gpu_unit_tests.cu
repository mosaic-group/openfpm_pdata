#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "Grid/grid_dist_id.hpp"

BOOST_AUTO_TEST_SUITE( sgrid_gpu_test_suite )

template<unsigned int p>
struct insert_kernel2D
{
	template<typename SparseGridGpu_type, typename ite_type>
	__device__ void operator()(SparseGridGpu_type & sg, ite_type & ite, float c)
	{
		GRID_ID_2_GLOBAL(ite);

	    sg.init();

	    sg.template insert<p>(key) = c + keyg.get(0) + keyg.get(1);

	    __syncthreads();

	    sg.flush_block_insert();
	}
};

template<unsigned int p>
struct insert_kernel3D
{
	template<typename SparseGridGpu_type, typename ite_type>
	__device__ void operator()(SparseGridGpu_type & sg, ite_type & ite, float c)
	{
		GRID_ID_3_GLOBAL(ite);

	    sg.init();

	    sg.template insert<p>(key) = c + keyg.get(0) + keyg.get(1) + keyg.get(2);

	    __syncthreads();

	    sg.flush_block_insert();
	}
};

template<unsigned int p>
struct stencil_kernel
{
	template<typename SparseGridGpu_type>
	__device__ void operator()(SparseGridGpu_type & sg, ite_gpu<SparseGridGpu_type::d> & ite, float c)
	{
		// TODO
	}
};

BOOST_AUTO_TEST_CASE( sgrid_gpu_test_base )
{
	size_t sz[2] = {17,17};
	periodicity<2> bc = {PERIODIC,PERIODIC};

	Ghost<2,long int> g(1);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	sgrid_dist_id_gpu<2,float,aggregate<float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{1,1});
	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	/////// GPU Run kernel

	float c = 5.0;

	gdist.template iterateGridGPU<insert_kernel2D<0>>(it,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template deviceToHost<0>();

	{
		Box<2,size_t> box2({0,0},{15,15});

		auto it = gdist.getGridIterator(box2.getKP1(),box2.getKP2());

		while (it.isNext())
		{
			auto p = it.get_dist();
			auto p2 = it.get();

			if (p2.get(0) == box.getLow(0) && p2.get(1) == box.getLow(1))
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 7.0);
			}
			else
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 666.0);
			}

			++it;
		}
	}

	//

	c = 3.0;

	Box<2,size_t> box3({3,3},{11,11});

	auto it3 = gdist.getGridIterator(box3.getKP1(),box3.getKP2());

	gdist.template iterateGridGPU<insert_kernel2D<0>>(it3,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);
	gdist.template deviceToHost<0>();

	{
		Box<2,size_t> box2({0,0},{15,15});

		auto it = gdist.getGridIterator(box2.getKP1(),box2.getKP2());

		while (it.isNext())
		{
			auto p = it.get_dist();
			auto p2 = it.get();

			Point<2,size_t> p2_ = p2.toPoint();

			if (box.isInside(p2_))
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 7.0);
			}
			else if (box3.isInside(p2_))
			{
				float tst = c + p2.get(0) + p2.get(1);

				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), tst);
			}
			else
			{
				BOOST_REQUIRE_EQUAL(gdist.template get<0>(p), 666.0);
			}

			++it;
		}
	}
}


BOOST_AUTO_TEST_CASE( sgrid_gpu_test_output )
{
	auto & v_cl = create_vcluster();

	if (v_cl.size() > 3){return;}

	size_t sz[2] = {17,17};
	periodicity<2> bc = {PERIODIC,PERIODIC};

	Ghost<2,long int> g(1);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	sgrid_dist_id_gpu<2,float,aggregate<float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{15,15});
	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	/////// GPU Run kernel

	float c = 5.0;

	gdist.template iterateGridGPU<insert_kernel2D<0>>(it,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template deviceToHost<0>();

	gdist.write("sgrid_gpu_output");

	std::string file_test("sgrid_gpu_output_" + std::to_string(v_cl.size()) + "_" + std::to_string(v_cl.rank())  + ".vtk");
	std::string file("sgrid_gpu_output_" + std::to_string(v_cl.rank()) + ".vtk");

	bool test = compare(file,"test_data/" + file_test);

	BOOST_REQUIRE_EQUAL(true,test);
}

void sgrid_ghost_get(size_t (& sz)[2],size_t (& sz2)[2])
{
	periodicity<2> bc = {PERIODIC,PERIODIC};

	Ghost<2,long int> g(1);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	sgrid_dist_id_gpu<2,float,aggregate<float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{sz2[0],sz2[1]});
	auto it = gdist.getGridIterator(box.getKP1(),box.getKP2());

	/////// GPU Run kernel

	float c = 5.0;

	gdist.template iterateGridGPU<insert_kernel2D<0>>(it,c);
	gdist.template flush<smax_<0>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template deviceToHost<0>();
	gdist.write_debug("before_ghost");

	gdist.template ghost_get<0>(RUN_ON_DEVICE);

	gdist.template deviceToHost<0>();
	gdist.write_debug("after_ghost");

	// Now we check that ghost is correct

	auto it2 = gdist.getDomainIterator();

	bool match = true;

	while (it2.isNext())
	{
		auto p = it2.get();

		auto key = it2.getGKey(p);

		auto p_xp1 = p.move(0,1);
		auto p_xm1 = p.move(0,-1);
		auto p_yp1 = p.move(1,1);
		auto p_ym1 = p.move(1,-1);

		auto key_xp1 = key.move(0,1);
		auto key_xm1 = key.move(0,-1);
		auto key_yp1 = key.move(1,1);
		auto key_ym1 = key.move(1,-1);

		if (box.isInside(key_xp1.toPoint()))
		{
			match &= gdist.template get<0>(p_xp1) == c + key_xp1.get(0) + key_xp1.get(1);

			if (match == false)
			{
				std::cout << gdist.template get<0>(p_xp1) << "   " << c + key_xp1.get(0) + key_xp1.get(1) << std::endl;
				break;
			}
		}

		if (box.isInside(key_xm1.toPoint()))
		{
			match &= gdist.template get<0>(p_xm1) == c + key_xm1.get(0) + key_xm1.get(1);

			if (match == false)
			{
				std::cout << gdist.template get<0>(p_xm1) << "   " << c + key_xm1.get(0) + key_xm1.get(1) << std::endl;
				break;
			}
		}

		if (box.isInside(key_yp1.toPoint()))
		{
			match &= gdist.template get<0>(p_yp1) == c + key_yp1.get(0) + key_yp1.get(1);

			if (match == false)
			{
				std::cout << gdist.template get<0>(p_yp1) << "   " << c + key_yp1.get(0) + key_yp1.get(1) << std::endl;
				break;
			}
		}

		if (box.isInside(key_ym1.toPoint()))
		{
			match &= gdist.template get<0>(p_ym1) == c + key_ym1.get(0) + key_ym1.get(1);

			if (match == false)
			{
				std::cout << gdist.template get<0>(p_ym1) << "   " << c + key_ym1.get(0) + key_ym1.get(1) << std::endl;
				break;
			}
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( sgrid_gpu_test_ghost_get )
{
	size_t sz[2] = {17,17};
	size_t sz6[2] = {15,15};
	sgrid_ghost_get(sz,sz6);

	size_t sz2[2] = {170,170};
	size_t sz3[2] = {15,15};
	sgrid_ghost_get(sz2,sz3);

	size_t sz4[2] = {168,168};
	sgrid_ghost_get(sz2,sz4);
}


BOOST_AUTO_TEST_CASE( sgrid_gpu_test_conv2_test )
{
	size_t sz[2] = {164,164};
	periodicity<2> bc = {PERIODIC,PERIODIC};

	Ghost<2,long int> g(1);

	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	sgrid_dist_id_gpu<2,float,aggregate<float,float,float,float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);
	gdist.template setBackgroundValue<1>(666);
	gdist.template setBackgroundValue<2>(666);
	gdist.template setBackgroundValue<3>(666);

	/////// GPU insert + flush

	Box<2,size_t> box({1,1},{sz[0],sz[1]});

	/////// GPU Run kernel

	float c = 5.0;

	typedef typename GetAddBlockType<decltype(gdist)>::type InsertBlockT;

	gdist.addPoints(box.getKP1(),box.getKP2(),
			        [] __device__ (int i, int j)
			        {
						return true;
			        },
			        [c] __device__ (InsertBlockT & data, int i, int j)
			        {
			        	data.template get<0>() = c + i + j;
			        	data.template get<1>() = c + 1000 + i + j;
			        }
			        );

	gdist.template flush<smax_<0>,smax_<1>>(flush_type::FLUSH_ON_DEVICE);
	gdist.template ghost_get<0,1>(RUN_ON_DEVICE);


	// Now run the convolution

	typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

	gdist.template conv2<0,1,2,3,1>({2,2},{(int)sz[0]-2,(int)sz[1]-2},[] __device__ (float & u_out, float & v_out, CpBlockType & u, CpBlockType & v,int i, int j){
		u_out = u(i+1,j) - u(i-1,j) + u(i,j+1) - u(i,j-1);
		v_out = v(i+1,j) - v(i-1,j) + v(i,j+1) - v(i,j-1);
	});

	gdist.deviceToHost<0,1,2,3>();

	// Now we check that ghost is correct

	auto it3 = gdist.getSubDomainIterator({2,2},{(int)sz[0]-2,(int)sz[1]-2});

	bool match = true;

	while (it3.isNext())
	{
		auto p = it3.get();

		auto p_xp1 = p.move(0,1);
		auto p_xm1 = p.move(0,-1);
		auto p_yp1 = p.move(1,1);
		auto p_ym1 = p.move(1,-1);

		float sub1 = gdist.template get<2>(p);
		float sub2 = gdist.template get<3>(p);

		if (sub1 != 4.0 || sub2 != 4.0)
		{
			std::cout << sub1 << "  " << sub2 << std::endl;
			std::cout << gdist.template get<0>(p_xp1) << "   " << gdist.template get<0>(p_xm1) << std::endl;
			std::cout << gdist.template get<1>(p_xp1) << "   " << gdist.template get<1>(p_xm1) << std::endl;
			match = false;
			break;
		}

		++it3;
	}

	gdist.write("SGRID");

	BOOST_REQUIRE_EQUAL(match,true);
}


BOOST_AUTO_TEST_CASE( sgrid_gpu_test_conv2_test_3d )
{
	size_t sz[3] = {60,60,60};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};

	Ghost<3,long int> g(1);

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	sgrid_dist_id_gpu<3,float,aggregate<float,float,float,float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);
	gdist.template setBackgroundValue<1>(666);
	gdist.template setBackgroundValue<2>(666);
	gdist.template setBackgroundValue<3>(666);

	/////// GPU insert + flush

	Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});

	/////// GPU Run kernel

	float c = 5.0;

	typedef typename GetAddBlockType<decltype(gdist)>::type InsertBlockT;

	gdist.addPoints(box.getKP1(),box.getKP2(),
			        [] __device__ (int i, int j, int k)
			        {
						return true;
			        },
			        [c] __device__ (InsertBlockT & data, int i, int j, int k)
			        {
			        	data.template get<0>() = c + i + j + k;
			        	data.template get<1>() = c + 1000 + i + j + k;
			        }
			        );

	gdist.template flush<smax_<0>,smax_<1>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template ghost_get<0,1>(RUN_ON_DEVICE);

	for (int i = 0 ; i < 10 ; i++)
	{
		gdist.template ghost_get<0,1>(RUN_ON_DEVICE);
	}

	// Now run the convolution

	typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

	gdist.template conv2<0,1,2,3,1>({2,2,2},{(int)sz[0]-2,(int)sz[1]-2,(int)sz[2]-2},[] __device__ (float & u_out, float & v_out, CpBlockType & u, CpBlockType & v,int i, int j, int k){
		u_out = u(i+1,j,k) - u(i-1,j,k) + u(i,j+1,k) - u(i,j-1,k) + u(i,j,k+1) - u(i,j,k-1);
		v_out = v(i+1,j,k) - v(i-1,j,k) + v(i,j+1,k) - v(i,j-1,k) + v(i,j,k+1) - v(i,j,k-1);
	});

	gdist.deviceToHost<0,1,2,3>();

	// Now we check that ghost is correct

	auto it3 = gdist.getSubDomainIterator({2,2,2},{(int)sz[0]-2,(int)sz[1]-2,(int)sz[2]-2});

	bool match = true;

	while (it3.isNext())
	{
		auto p = it3.get();

		auto p_xp1 = p.move(0,1);
		auto p_xm1 = p.move(0,-1);
		auto p_yp1 = p.move(1,1);
		auto p_ym1 = p.move(1,-1);
		auto p_zp1 = p.move(2,1);
		auto p_zm1 = p.move(2,-1);

		float sub1 = gdist.template get<2>(p);
		float sub2 = gdist.template get<3>(p);

		if (sub1 != 6.0 || sub2 != 6.0)
		{
			std::cout << sub1 << "  " << sub2 << std::endl;
			std::cout << gdist.template get<0>(p_xp1) << "   " << gdist.template get<0>(p_xm1) << std::endl;
			std::cout << gdist.template get<1>(p_xp1) << "   " << gdist.template get<1>(p_xm1) << std::endl;
			match = false;
			break;
		}

		++it3;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( sgrid_gpu_test_ghost_point_remove )
{
	size_t sz[3] = {60,60,60};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};

	Ghost<3,long int> g(1);

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	sgrid_dist_id_gpu<3,float,aggregate<float,float,float,float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);
	gdist.template setBackgroundValue<1>(666);
	gdist.template setBackgroundValue<2>(666);
	gdist.template setBackgroundValue<3>(666);

	/////// GPU insert + flush

	Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});

	/////// GPU Run kernel

	float c = 5.0;

	typedef typename GetAddBlockType<decltype(gdist)>::type InsertBlockT;

	gdist.addPoints(box.getKP1(),box.getKP2(),
			        [] __device__ (int i, int j, int k)
			        {
						return true;
			        },
			        [c] __device__ (InsertBlockT & data, int i, int j, int k)
			        {
			        	data.template get<0>() = c + i + j + k;
			        	data.template get<1>() = c + 1000 + i + j + k;
			        }
			        );

	gdist.template flush<smax_<0>,smax_<1>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template ghost_get<0,1>(RUN_ON_DEVICE);

	// Remove the right side of the points
	Box<3,size_t> bxR({59,0,0},{59,59,59});
	gdist.removePoints(bxR);

	// Remove the right side of the points
	Box<3,size_t> bxT({0,0,59},{59,59,59});
	gdist.removePoints(bxT);

	// Remove the right side of the points
	Box<3,size_t> bxD({0,59,0},{59,59,59});
	gdist.removePoints(bxD);

	gdist.template ghost_get<0,1>(RUN_ON_DEVICE);

	for (int i = 0 ; i < 10 ; i++)
	{
		gdist.template ghost_get<0,1>(RUN_ON_DEVICE);
	}

	// Now run the convolution

	typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

	gdist.template conv2<0,1,2,3,1>({2,2,2},{(int)sz[0]-3,(int)sz[1]-3,(int)sz[2]-3},[] __device__ (float & u_out, float & v_out, CpBlockType & u, CpBlockType & v,int i, int j, int k){
		u_out = u(i+1,j,k) - u(i-1,j,k) + u(i,j+1,k) - u(i,j-1,k) + u(i,j,k+1) - u(i,j,k-1);
		v_out = v(i+1,j,k) - v(i-1,j,k) + v(i,j+1,k) - v(i,j-1,k) + v(i,j,k+1) - v(i,j,k-1);
	});

	gdist.deviceToHost<0,1,2,3>();

	// Now we check that ghost is correct

	auto it3 = gdist.getSubDomainIterator({2,2,2},{(int)sz[0]-3,(int)sz[1]-3,(int)sz[2]-3});

	bool match = true;

	while (it3.isNext())
	{
		auto p = it3.get();

		auto p_xp1 = p.move(0,1);
		auto p_xm1 = p.move(0,-1);
		auto p_yp1 = p.move(1,1);
		auto p_ym1 = p.move(1,-1);
		auto p_zp1 = p.move(2,1);
		auto p_zm1 = p.move(2,-1);

		float sub1 = gdist.template get<2>(p);
		float sub2 = gdist.template get<3>(p);

		if (sub1 != 6.0 || sub2 != 6.0)
		{
			std::cout << sub1 << "  " << sub2 << std::endl;
			std::cout << gdist.template get<0>(p_xp1) << "   " << gdist.template get<0>(p_xm1) << std::endl;
			std::cout << gdist.template get<1>(p_xp1) << "   " << gdist.template get<1>(p_xm1) << std::endl;
			match = false;
			break;
		}

		++it3;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	gdist.template deviceToHost<0,1,2,3>();

	auto it4 = gdist.getDomainGhostIterator();
	Box<3,long int> bin({0,0,0},{59,59,59});

	match = true;

	while (it4.isNext())
	{
		auto p = it4.get();

		// We have to check we have no point in the ghost area
		auto gkey = it4.getGKey(p);

		if (bin.isInside(gkey.toPoint()) == false)
		{match = false;}

		++it4;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( sgrid_gpu_test_skip_labelling )
{
	size_t sz[3] = {60,60,60};
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};

	Ghost<3,long int> g(1);

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	sgrid_dist_id_gpu<3,float,aggregate<float,float,float,float>> gdist(sz,domain,g,bc);

	gdist.template setBackgroundValue<0>(666);
	gdist.template setBackgroundValue<1>(666);
	gdist.template setBackgroundValue<2>(666);
	gdist.template setBackgroundValue<3>(666);

	/////// GPU insert + flush

	Box<3,size_t> box({1,1,1},{sz[0]-1,sz[1]-1,sz[2]-1});

	/////// GPU Run kernel

	float c = 5.0;

	typedef typename GetAddBlockType<decltype(gdist)>::type InsertBlockT;

	gdist.addPoints(box.getKP1(),box.getKP2(),
			        [] __device__ (int i, int j, int k)
			        {
						return true;
			        },
			        [c] __device__ (InsertBlockT & data, int i, int j, int k)
			        {
			        	data.template get<0>() = c + i + j + k;
			        	data.template get<1>() = c + 1000 + i + j + k;
			        }
			        );

	gdist.template flush<smax_<0>,smax_<1>>(flush_type::FLUSH_ON_DEVICE);

	gdist.template ghost_get<0,1>(RUN_ON_DEVICE);

	// Now run the convolution

	typedef typename GetCpBlockType<decltype(gdist),0,1>::type CpBlockType;

	gdist.template conv2<0,1,0,1,1>({0,0,0},{(int)sz[0]-1,(int)sz[1]-1,(int)sz[2]-1},[] __device__ (float & u_out, float & v_out, CpBlockType & u, CpBlockType & v,int i, int j, int k){
		u_out = 1*u(i,j,k);
		v_out = 1*v(i,j,k);
	});

	gdist.template ghost_get<0,1>(RUN_ON_DEVICE | SKIP_LABELLING);

	gdist.template conv2<0,1,0,1,1>({0,0,0},{(int)sz[0]-1,(int)sz[1]-1,(int)sz[2]-1},[] __device__ (float & u_out, float & v_out, CpBlockType & u, CpBlockType & v,int i, int j, int k){
		u_out = 5*u(i,j,k);
		v_out = 5*v(i,j,k);
	});

	gdist.template ghost_get<0,1>(RUN_ON_DEVICE | SKIP_LABELLING);

	//////////////////////////////////// DEBUG ///////////////////////////////
	gdist.deviceToHost<0,1,2,3>();
	gdist.write_debug("DEBUG");
	//////////////////////////////////////////////////////////////////////////

	gdist.template conv2<0,1,0,1,1>({0,0,0},{(int)sz[0]-1,(int)sz[1]-1,(int)sz[2]-1},[] __device__ (float & u_out, float & v_out, CpBlockType & u, CpBlockType & v,int i, int j, int k){
		u_out = 2*u(i,j,k);
		v_out = 2*v(i,j,k);
	});

	gdist.template ghost_get<0,1>(RUN_ON_DEVICE | SKIP_LABELLING);

	gdist.template conv2<0,1,2,3,1>({2,2,2},{(int)sz[0]-3,(int)sz[1]-3,(int)sz[2]-3},[] __device__ (float & u_out, float & v_out, CpBlockType & u, CpBlockType & v,int i, int j, int k){
		u_out = u(i+1,j,k) - u(i-1,j,k) + u(i,j+1,k) - u(i,j-1,k) + u(i,j,k+1) - u(i,j,k-1);
		v_out = v(i+1,j,k) - v(i-1,j,k) + v(i,j+1,k) - v(i,j-1,k) + v(i,j,k+1) - v(i,j,k-1);
	});


	gdist.deviceToHost<0,1,2,3>();

	// Now we check that ghost is correct

	auto it3 = gdist.getSubDomainIterator({2,2,2},{(int)sz[0]-3,(int)sz[1]-3,(int)sz[2]-3});

	bool match = true;

	while (it3.isNext())
	{
		auto p = it3.get();

		auto p_xp1 = p.move(0,1);
		auto p_xm1 = p.move(0,-1);
		auto p_yp1 = p.move(1,1);
		auto p_ym1 = p.move(1,-1);
		auto p_zp1 = p.move(2,1);
		auto p_zm1 = p.move(2,-1);

		float sub1 = gdist.template get<2>(p);
		float sub2 = gdist.template get<3>(p);

		if (sub1 != 6.0*10.0 || sub2 != 6.0*10.0)
		{
			std::cout << sub1 << "  " << sub2 << std::endl;
			std::cout << gdist.template get<0>(p_xp1) << "   " << gdist.template get<0>(p_xm1) << std::endl;
			std::cout << gdist.template get<1>(p_xp1) << "   " << gdist.template get<1>(p_xm1) << std::endl;
			match = false;
			break;
		}

		++it3;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_SUITE_END()
