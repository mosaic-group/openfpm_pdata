 /*! \page Vector Vector
 *
 * \subpage Vector_0_simple
 * \subpage Vector_1_celllist
 * \subpage Vector_1_ghost_get
 * \subpage Vector_1_HDF5
 * \subpage Vector_1_gpu_first_step
 * \subpage Vector_2_expression
 * \subpage Vector_3_md
 * \subpage Vector_3_md_dyn_gpu
 * \subpage Vector_3_md_dyn_gpu_opt
 * \subpage Vector_4_reo_root
 * \subpage Vector_4_cp
 * \subpage Vector_4_mp_cl
 * \subpage Vector_5_md_vl_sym
 * \subpage Vector_5_md_vl_sym_crs
 * \subpage Vector_6_complex_usage
 * \subpage Vector_7_sph_dlb
 * \subpage Vector_7_sph_dlb_opt
 * \subpage Vector_7_sph_dlb_gpu
 * \subpage Vector_7_sph_dlb_gpu_opt
 * \subpage Vector_8_DEM
 * \subpage Vector_9_gpu_cuda_interop
 *
 */


/*!
 * \page Vector_0_simple Vector 0 simple
 *
 *
 * [TOC]
 *
 *
 * # Simple Vector example # {#simple_vector_example}
 *
 *
 * This example show several basic functionalities of the distributed vector, A distributed vector is nothing else than
 * a set of particles in an N dimensional space
 *
 * \htmlonly
 * <a href="#" onclick="hide_show('vector-video-1')" >Video 1</a>
 * <div style="display:none" id="vector-video-1">
 * <video id="vid1" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-1.mp4" type="video/mp4"></video>
 * <script>video_anim('vid1',100,230)</script>
 * </div>
 * <a href="#" onclick="hide_show('vector-video-2')" >Video 2</a>
 * <div style="display:none" id="vector-video-2">
 * <video id="vid2" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-2.mp4" type="video/mp4"></video>
 * <script>video_anim('vid2',21,1590)</script>
 * </div>
 * \endhtmlonly
 *
 * ## inclusion ## {#e0_v_inclusion}
 *
 * In order to use distributed vectors in our code we have to include the file Vector/vector_dist.hpp
 *
 * \snippet Vector/0_simple/main.cpp inclusion
 *
 */

//! \cond [inclusion] \endcond
#include <stddef.h>
#include "Vector/vector_dist.hpp"
//! \cond [inclusion] \endcond

int main(int argc, char* argv[])
{

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Initialization ## {#e0_s_init}
	 *
	 *  Here we
	 *  * Initialize the library
	 *  * we create a Box that define our domain
	 *  * An array that define our boundary conditions
	 *  * A Ghost object that will define the extension of the ghost part in physical units
	 *
	 *
	 *
	 * \snippet Vector/0_simple/main.cpp Initialization and parameters
	 *
	 * \htmlonly
	 * <a href="#" onclick="hide_show('vector-video-5')" >Video 1</a>
	 * <div style="display:none" id="vector-video-5">
	 * <video id="vid5" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-4.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid5',447,513)</script>
	 * </div>
	 * <a href="#" onclick="hide_show('vector-video-4')" >Video 2</a>
	 * <div style="display:none" id="vector-video-4">
	 * <video id="vid4" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-4.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid4',594,1023)</script>
	 * </div>
	 * \endhtmlonly
	 *
	 */

	//! \cond [Initialization and parameters] \endcond

    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.01);
	
	//! \cond [Initialization and parameters] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Vector instantiation ## {#e0_s_vector_inst}
	 *
	 * Here we are creating a distributed vector defined by the following parameters
	 *
	 * * 2 is the Dimensionality of the space where the objects live
	 * * float is the type used for the spatial coordinate of the particles
	 * * float,float[3],float[3][3] is the information stored by each particle a scalar float, a vector float[3] and a tensor of rank 2 float[3][3]
	 *   the list of properties must be put into an aggregate data structure aggregate<prop1,prop2,prop3, ... >
	 *
	 * vd is the instantiation of the object
	 *
	 * The Constructor instead require:
	 *
	 * * Number of particles 4096 in this case
	 * * Domain where is defined this structure
	 * * bc boundary conditions
	 * * g Ghost
	 *
	 * The following construct a vector where each processor has 4096 / N_proc (N_proc = number of processor)
	 * objects with an undefined position in space. This non-space decomposition is also called data-driven
	 * decomposition
	 *
	 *
	 * \snippet Vector/0_simple/main.cpp vector instantiation
	 *
	 * \htmlonly
	 * <a href="#" onclick="hide_show('vector-video-3')" >Video</a>
	 * <div style="display:none" id="vector-video-3">
	 * <video id="vid3" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-4.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid3',1047,1370)</script>
	 * </div>
	 * \endhtmlonly
	 *
	 */

	//! \cond [vector instantiation] \endcond

    vector_dist<2,float, aggregate<float,float[3],float[3][3]> > vd(4096,domain,bc,g);

	// the scalar is the element at position 0 in the aggregate
	const int scalar = 0;

	// the vector is the element at position 1 in the aggregate
	const int vector = 1;

	// the tensor is the element at position 2 in the aggregate
	const int tensor = 2;

	//! \cond [vector instantiation] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Assign position ## {#e0_s_assign_pos}
	 *
	 * Get an iterator that go through the 4096 particles. Initially all the particles
	 *  has an undefined position state. In this cycle we define its position. In this
	 * example we use iterators. Iterators are convenient way to explore/iterate data-structures in an
	 * convenient and easy way
	 *
	 *  \snippet Vector/0_simple/main.cpp assign position
	 *
	 * \htmlonly
	 * <a href="#" onclick="hide_show('vector-video-13')" >Iterator Video</a>
	 * <div style="display:none" id="vector-video-13">
	 * <video id="vid13" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-8.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid13',31,1362)</script>
	 * </div>
	 * \endhtmlonly
	 *
	 */

	//! \cond [assign position] \endcond

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(key)[0] = (float)rand() / RAND_MAX;

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(key)[1] = (float)rand() / RAND_MAX;

		// next particle
		++it;
	}

	//! \cond [assign position] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Mapping particles ## {#e0_s_map}
	 *
	 * On a parallel program, once we define the position, we distribute the particles according to the underlying space decomposition
	 * The default decomposition is created even before assigning the position to the object, and is calculated
	 * giving to each processor an equal portion of space minimizing the surface to reduce communication.
	 *
	 * \snippet Vector/0_simple/main.cpp map
	 *
	 * \htmlonly
	 * <a href="#" onclick="hide_show('vector-video-11')" >Parallelization Video</a>
	 * <div style="display:none" id="vector-video-11">
	 * <video id="vid11" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-5.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid11',440,995)</script>
	 * </div>
	 * <a href="#" onclick="hide_show('vector-video-8')" >Video 1</a>
	 * <div style="display:none" id="vector-video-8">
	 * <video id="vid8" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-6.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid8',0,483)</script>
	 * </div>
	 * <a href="#" onclick="hide_show('vector-video-9')" >Video 2</a><br>
	 * <div style="display:none" id="vector-video-9">
	 * <video id="vid9" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-6.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid9',1009,1041)</script>
	 * </div>
	 * <a href="#" onclick="hide_show('vector-video-10')" >Video 3</a>
	 * <div style="display:none" id="vector-video-10">
	 * <video id="vid10" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-6.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid10',1739,1948)</script>
	 * </div>
	 * \endhtmlonly
	 *
	 */

	//! \cond [map] \endcond

	vd.map();

	//! \cond [map] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Assign values to particles property ## {#assign_prop}
	 *
	 * We Iterate across all the particles, we count them using a local counter and we assign 1.0 to
	 * all the particles properties. Each particle has a scalar, vector and tensor property.
	 *
	 * \snippet Vector/0_simple/main.cpp assign property
	 *
	 *
	 */

	//! \cond [assign property] \endcond

	//Counter we use it later
	size_t cnt = 0;

	// Get a particle iterator
	it = vd.getDomainIterator();

	// For each particle ...
	while (it.isNext())
	{
		// ... p
		auto p = it.get();

		// we set the properties of the particle p
		
         // the scalar property
		vd.template getProp<scalar>(p) = 1.0;

		vd.template getProp<vector>(p)[0] = 1.0;
		vd.template getProp<vector>(p)[1] = 1.0;
		vd.template getProp<vector>(p)[2] = 1.0;

		vd.template getProp<tensor>(p)[0][0] = 1.0;
		vd.template getProp<tensor>(p)[0][1] = 1.0;
		vd.template getProp<tensor>(p)[0][2] = 1.0;
		vd.template getProp<tensor>(p)[1][0] = 1.0;
		vd.template getProp<tensor>(p)[1][1] = 1.0;
		vd.template getProp<tensor>(p)[1][2] = 1.0;
		vd.template getProp<tensor>(p)[2][0] = 1.0;
		vd.template getProp<tensor>(p)[2][1] = 1.0;
        vd.template getProp<tensor>(p)[2][2] = 1.0;

		// increment the counter
		cnt++;

		// next particle
		++it;
	}

	//! \cond [assign property] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Reduce (sum numbers across processors) ## {#e0_s_reduce}
	 *
	 * cnt contain the number of object the local processor contain, if we are interested to count the total number across the processors
	 * we can use the function add, to sum across the processors. First we have to get an instance of Vcluster, queue an operation of add with
	 * the variable count and finally execute. All the operations are asynchronous, execute work like a barrier and ensure that all the
	 * queued operations are executed.
	 *
	 * \snippet Vector/0_simple/main.cpp reduce
	 *
	 */

	//! \cond [reduce] \endcond
	
	auto & v_cl = create_vcluster();
	v_cl.sum(cnt);
	v_cl.execute();
	
	//! \cond [reduce] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Visualization, write VTK files ## {#e0_s_vis_vtk}
	 *
	 * With this function we output the particle position in VTK format. A VTK file
	 * contain information about particle position and properties. Such file can be visualized
	 * with a program like paraview. In case this program run on several processor each processor
	 * generate a VTK file. VTK has two format one is the ASCII that is human readable but produce
	 * bigger file the other is the binary that produce not-human readable files, but smaller
	 * and more efficent files to read for Paraview. In order to create a Binary VTK file use the
	 * option VTK_WRITER in combination with FORMAT_BINARY. By default properties does not have name
	 * the VTK writer will label them "attr0" for the first property "attr1" for the second and so on
	 * in case we want to assign names to the properties we can use the function setPropNames
	 *
	 * \snippet Vector/0_simple/main.cpp vtk
	 *
	 * \htmlonly
	 * <a href="#" onclick="hide_show('vector-video-6')" >Video</a>
	 * <div style="display:none" id="vector-video-6">
	 * <video id="vid6" width="1200" height="576" controls> <source src="http://ppmcore.mpi-cbg.de/upload/video/Lesson1-5.mp4" type="video/mp4"></video>
	 * <script>video_anim('vid6',92,400)</script>
	 * </div>
	 * \endhtmlonly
	 *
	 */

	//! \cond [vtk] \endcond

	openfpm::vector<std::string> names({"scalar","vector","tensor"});
	vd.setPropNames(names);

	// save vtk format (vtk is always the default)
	vd.write("particles");

	// save in vtk binary format
	vd.write("particles_bin",VTK_WRITER | FORMAT_BINARY);

	//! \cond [vtk] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Finalize ## {#finalize_e0_sim}
	 *
	 *  At the very end of the program we have always de-initialize the library
	 *
	 * \snippet Vector/0_simple/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Full code ## {#code_e0_sim}
	 *
	 * \include Vector/0_simple/main.cpp
	 *
	 */
}
