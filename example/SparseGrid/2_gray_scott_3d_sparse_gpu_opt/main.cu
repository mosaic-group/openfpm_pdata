//#define VCLUSTER_PERF_REPORT <------ Activate telemetry for the VCluster data-structure
//#define SYNC_BEFORE_TAKE_TIME <------ Force synchronization of the kernels everytime we take the time with the structure timer.
//                                      Use this option for telemetry and GPU otherwise the result are unreliable                                        
//#define ENABLE_GRID_DIST_ID_PERF_STATS  <------ Activate telementry for the grid data-structure

#include "Decomposition/Distribution/BoxDistribution.hpp"
#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "timer.hpp"

/*!
 *
 * \page Grid_3_gs_3D_sparse_gpu_opt Gray Scott in 3D using sparse grids on GPU (Optimized)
 *
 * [TOC]
 *
 * # Solving a gray scott-system in 3D using Sparse grids on gpu (Optimized) # {#e3_gs_gray_scott_gpu}
 *
 * This example show how to solve a Gray-Scott system in 3D using sparse grids on gpu
 *
 * In figure is the final solution of the problem
 *
 * \htmlonly
 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/gray_scott_3d/gs_alpha.png"/>
 * \endhtmlonly
 *
 * More or less this example is the adaptation of the dense example in 3D
 *
 * \see \ref Grid_3_gs_3D
 *
 * # Optimizations
 *
 * Instead of using the default decomposition algorithm based on parmetis we use BoxDistribution. This decomposition divide the space equally 
 * across processors. The way to use a different algorithm for decomposing the sparse grid is given by changing the type of the Sparse grid
 * 
 * \snippet SparseGrid/2_gray_scott_3d_sparse_gpu_opt/main.cu grid definition
 *
 * Because the geometry is fixed we are also using the option SKIP_LABELLING. With this option active after a normal ghost_get we are able to 
 * activate certain optimization patterns in constructions of the sending buffers and merging data.
 *
 */

#ifdef __NVCC__

constexpr int U = 0;
constexpr int V = 1;

constexpr int U_next = 2;
constexpr int V_next = 3;

constexpr int x = 0;
constexpr int y = 1;
constexpr int z = 2;

//! \cond [grid definition] \endcond

typedef CartDecomposition<3,float, CudaMemory, memory_traits_inte, BoxDistribution<3,float> > Dec;

typedef sgrid_dist_id_gpu<3,float,aggregate<float,float,float,float>,CudaMemory, Dec> SparseGridType;

//! \cond [grid definition] \endcond

void init(SparseGridType & grid, Box<3,float> & domain)
{
	//! \cond [create points] \endcond

	typedef typename GetAddBlockType<SparseGridType>::type InsertBlockT;

	grid.addPoints([] __device__ (int i, int j, int k)
			        {
						return true;
			        },
			        [] __device__ (InsertBlockT & data, int i, int j, int k)
			        {
			        	data.template get<U>() = 1.0;
			        	data.template get<V>() = 0.0;
			        }
			        );


	grid.template flush<smax_<U>,smax_<V>>(flush_type::FLUSH_ON_DEVICE);

	//! \cond [create points] \endcond

	long int x_start = grid.size(0)*1.55f/domain.getHigh(0);
	long int y_start = grid.size(1)*1.55f/domain.getHigh(1);
	long int z_start = grid.size(1)*1.55f/domain.getHigh(2);

	long int x_stop = grid.size(0)*1.85f/domain.getHigh(0);
	long int y_stop = grid.size(1)*1.85f/domain.getHigh(1);
	long int z_stop = grid.size(1)*1.85f/domain.getHigh(2);

	//! \cond [create points sub] \endcond

	grid_key_dx<3> start({x_start,y_start,z_start});
	grid_key_dx<3> stop ({x_stop,y_stop,z_stop});

        grid.addPoints(start,stop,[] __device__ (int i, int j, int k)
                                {
                                                return true;
                                },
                                [] __device__ (InsertBlockT & data, int i, int j, int k)
                                {
                                        data.template get<U>() = 0.5;
                                        data.template get<V>() = 0.24;
                                }
                                );

	grid.template flush<smax_<U>,smax_<V>>(flush_type::FLUSH_ON_DEVICE);

	//! \cond [create points sub] \endcond
}


int main(int argc, char* argv[])
{
	openfpm_init(&argc,&argv);

	// domain
	Box<3,float> domain({0.0,0.0,0.0},{2.5,2.5,2.5});
	
	// grid size
        size_t sz[3] = {256,256,256};

	// Define periodicity of the grid
	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};
	
	// Ghost in grid unit
	Ghost<3,long int> g(1);
	
	// deltaT
	float deltaT = 0.25;

	// Diffusion constant for specie U
	float du = 2*1e-5;

	// Diffusion constant for specie V
	float dv = 1*1e-5;

	// Number of timesteps
#ifdef TEST_RUN
	size_t timeSteps = 300;
#else
        size_t timeSteps = 15000;
#endif

	// K and F (Physical constant in the equation)
    float K = 0.053;
    float F = 0.014;

	SparseGridType grid(sz,domain,g,bc);

	// spacing of the grid on x and y
	float spacing[3] = {grid.spacing(0),grid.spacing(1),grid.spacing(2)};

	init(grid,domain);

	// sync the ghost
	grid.template ghost_get<U,V>(RUN_ON_DEVICE);

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2
	float uFactor = deltaT * du/(spacing[x]*spacing[x]);
	float vFactor = deltaT * dv/(spacing[x]*spacing[x]);

	auto & v_cl = create_vcluster();

	timer tot_sim;
	tot_sim.start();

	for (size_t i = 0; i < timeSteps ; ++i)
	{
		if (v_cl.rank() == 0)
		{std::cout << "STEP: " << i << std::endl;}
/*		if (i % 300 == 0)
		{
			std::cout << "STEP: " << i << std::endl;
			grid.write_frame("out",i,VTK_WRITER);
		}*/

		//! \cond [stencil get and use] \endcond

		typedef typename GetCpBlockType<decltype(grid),0,1>::type CpBlockType;

		//! \cond [lambda] \endcond

		auto func = [uFactor,vFactor,deltaT,F,K] __device__ (float & u_out, float & v_out,
				                                   CpBlockType & u, CpBlockType & v,
				                                   int i, int j, int k){

				float uc = u(i,j,k);
				float vc = v(i,j,k);

				u_out = uc + uFactor *(u(i-1,j,k) + u(i+1,j,k) +
                                                       u(i,j-1,k) + u(i,j+1,k) +
                                                       u(i,j,k-1) + u(i,j,k+1) - 6.0f*uc) - deltaT * uc*vc*vc
                                                       - deltaT * F * (uc - 1.0f);


				v_out = vc + vFactor *(v(i-1,j,k) + v(i+1,j,k) +
                                                       v(i,j+1,k) + v(i,j-1,k) +
                                                       v(i,j,k-1) + v(i,j,k+1) - 6.0f*vc) + deltaT * uc*vc*vc
					               - deltaT * (F+K) * vc;
				};

		//! \cond [lambda] \endcond

		//! \cond [body] \endcond

		if (i % 2 == 0)
		{
			cudaDeviceSynchronize();
			timer tconv;
			tconv.start();
			grid.conv2<U,V,U_next,V_next,1>({0,0,0},{(long int)sz[0]-1,(long int)sz[1]-1,(long int)sz[2]-1},func);
			cudaDeviceSynchronize();
			tconv.stop();
			std::cout << "Conv " << tconv.getwct() << std::endl;

			// After copy we synchronize again the ghost part U and V

			grid.ghost_get<U_next,V_next>(RUN_ON_DEVICE | SKIP_LABELLING);
		}
		else
		{
			grid.conv2<U_next,V_next,U,V,1>({0,0,0},{(long int)sz[0]-1,(long int)sz[1]-1,(long int)sz[2]-1},func);

			// After copy we synchronize again the ghost part U and V
			grid.ghost_get<U,V>(RUN_ON_DEVICE | SKIP_LABELLING);
		}

		//! \cond [body] \endcond

		// Every 500 time step we output the configuration for
		// visualization
//		if (i % 500 == 0)
//		{
//			grid.save("output_" + std::to_string(count));
//			count++;
//		}
	}
	
	tot_sim.stop();
	std::cout << "Total simulation: " << tot_sim.getwct() << std::endl;

	grid.deviceToHost<U,V,U_next,V_next>();
	grid.write("final");

	//! \cond [time stepping] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D
	 *
	 * ## Finalize ##
	 *
	 * Deinitialize the library
	 *
	 * \snippet Grid/3_gray_scott/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Grid_3_gs_3D_sparse Gray Scott in 3D
	 *
	 * # Full code # {#code}
	 *
	 * \include Grid/3_gray_scott_3d/main.cpp
	 *
	 */
}

#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif

