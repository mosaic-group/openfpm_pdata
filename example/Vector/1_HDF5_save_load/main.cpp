/*!
 * \page Vector_1_HDF5 Vector 1 HDF5 save and load
 *
 *
 * [TOC]
 *
 *
 * # HDF5 Save and load # {#HDF5_save_and_load_parallel}
 *
 *
 * This example show how to save and load a vector in the parallel
 *  format HDF5.
 *
 * ## inclusion ## {#e0_v_inclusion}
 *
 * In order to use distributed vectors in our code we have to include the file Vector/vector_dist.hpp
 *
 * \snippet Vector/1_HDF5_save_load/main.cpp inclusion
 *
 */

//! \cond [inclusion] \endcond

#include "Vector/vector_dist.hpp"

//! \cond [inclusion] \endcond

int main(int argc, char* argv[])
{

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
	 *
	 * ## Initialization ## {#HDF5_initialization}
	 *
	 *  Here we
	 *  * Initialize the library
	 *  * we create a Box that define our domain
	 *  * An array that define our boundary conditions
	 *  * A Ghost object that will define the extension of the ghost part in physical units
	 *
	 * \snippet Vector/1_HDF5_save_load/main.cpp Initialization and parameters
	 *
	 *
	 */

	//! \cond [Initialization and parameters] \endcond

    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<3,float> domain({0.0,0.0,0.0},{22.0,5.0,5.0});

	// Here we define the boundary conditions of our problem
    size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<3,float> g(0.01);
	
	//! \cond [Initialization and parameters] \endcond

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
	 *
	 * ## Vector instantiation ## {#HDF5_vector_inst}
	 *
	 * Here we are creating a distributed vector defined by the following parameters
	 *
	 * * **3** is the Dimensionality of the space where the objects live
	 * * **float** is the type used for the spatial coordinate of the particles
	 * * the information stored by each particle **float[3],float[3],float[3],float[3],float[3]**.
	 *   The list of properties must be put into an aggregate data structure aggregate<prop1,prop2,prop3, ... >
	 *
	 * vd is the instantiation of the object
	 *
	 * The Constructor instead require:
	 *
	 * * Number of particles, 4096 in this case
	 * * Domain where is defined this structure
	 * * bc boundary conditions
	 * * g Ghost
	 *
	 *
	 * \snippet Vector/1_HDF5_save_load/main.cpp vector instantiation
	 *
	 */

	//! \cond [vector instantiation] \endcond

	vector_dist<3,float, aggregate<float[3],float[3],float[3],float[3],float[3]> > vd(4096,domain,bc,g);

	//! \cond [vector instantiation] \endcond

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
	 *
	 * ## Assign position ## {#HDF5_vector_assign}
	 *
	 * Get an iterator that go through the 4096 particles. Initially all the particles
	 *  has an undefined position state. In this cycle we define its position. In this
	 * example we use iterators. Iterators are convenient way to explore/iterate data-structures in an
	 * convenient and easy way
	 *
	 * \snippet Vector/1_HDF5_save_load/main.cpp assign position
	 *
	 *
	 */

	//! \cond [assign position] \endcond

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(key)[0] = 22.0*((float)rand() / RAND_MAX);

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(key)[1] = 5.0*((float)rand() / RAND_MAX);

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(key)[1] = 5.0*((float)rand() / RAND_MAX);

		// next particle
		++it;
	}

	//! \cond [assign position] \endcond

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
	 *
	 * ## Mapping particles ## {#e0_s_map}
	 *
	 * On a parallel program, once we define the position, we distribute the particles according to the underlying space decomposition
	 * The default decomposition is created even before assigning the position to the object, and is calculated
	 * giving to each processor an equal portion of space minimizing the surface to reduce communication.
	 *
	 * \snippet Vector/1_HDF5_save_load/main.cpp map
	 *
	 *
	 */

	//! \cond [map] \endcond

	vd.map();

	//! \cond [map] \endcond

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
	 *
	 * ## Assign values to particles property ## {#assign_prop}
	 *
	 * We Iterate across all the particles, we assign some value
	 * to all the particles properties. Each particle has a scalar,
	 *  vector and tensor property.
	 *
	 * \snippet Vector/1_HDF5_save_load/main.cpp assign property
	 *
	 *
	 */

	//! \cond [assign property] \endcond

	// Get a particle iterator
	it = vd.getDomainIterator();

	// For each particle ...
	while (it.isNext())
	{
		// ... p
		auto p = it.get();

		// we set the properties of the particle p

		vd.template getProp<0>(p)[0] = 1.0;
		vd.template getProp<0>(p)[1] = 1.0;
		vd.template getProp<0>(p)[2] = 1.0;

		vd.template getProp<1>(p)[0] = 2.0;
		vd.template getProp<1>(p)[1] = 2.0;
		vd.template getProp<1>(p)[2] = 2.0;

		vd.template getProp<2>(p)[0] = 3.0;
		vd.template getProp<2>(p)[1] = 3.0;
		vd.template getProp<2>(p)[2] = 3.0;

		vd.template getProp<3>(p)[0] = 4.0;
		vd.template getProp<3>(p)[1] = 4.0;
		vd.template getProp<3>(p)[2] = 4.0;

		vd.template getProp<4>(p)[0] = 5.0;
		vd.template getProp<4>(p)[1] = 5.0;
		vd.template getProp<4>(p)[2] = 5.0;


		// next particle
		++it;
	}

	//! \cond [assign property] \endcond

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
	 *
	 * ## Parallel IO save ## {#HDF5_par_save}
	 *
	 * To save the file we use the function **save**. The system save
	 * all the information about the particles in the file (positions and
	 * all the properties, even complex properties)
	 *
	 * \see \ref vector_example_cp
	 *
	 * It is good to notice that independently from the number of
	 * processor this function produce only one file.
	 *
	 *
	 * \note The saved file can be used for checkpoint-restart. the status
	 *       of the particles or the application in general can be saved periodically
	 *        and can be used later-on to restart from the latest save
	 *
	 * \snippet Vector/1_HDF5_save_load/main.cpp save_part
	 *
	 */

	//! \cond [save_part] \endcond
	
	vd.save("particles_save.hdf5");

	//! \cond [save_part] \endcond

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
	 *
	 * ## Parallel IO load ## {#HDF5_par_load}
	 *
	 * To load the file we use the function **load**. The system load
	 * all the information about the particles (position and
	 * all the properties, even complex)
	 *
	 * \see \ref vector_example_cp
	 *
	 * It is good to notice that this function work independently
	 * from the number of processors used to save the file
	 *
	 * \warning Despite the fact that the function works for any number of processors
	 *          the **domain** parameter, the **dimensionality**, the **space type**
	 *          and the **properties** of the particles must match the saved one.
	 *          At the moment there is not check or control, so in case the
	 *          parameters does not match the load will produce error or
	 *          corrupted data.
	 *
	 *
	 * \snippet Vector/1_HDF5_save_load/main.cpp hdf5_load
	 *
	 */

	//! \cond [hdf5_load] \endcond

	vector_dist<3,float, aggregate<float[3],float[3],float[3],float[3],float[3]> > vd2(vd.getDecomposition(),0);

	// load the particles on another vector
	vd2.load("particles_save.hdf5");

	//! \cond [hdf5_load] \endcond

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
	 *
	 * ## Process loaded data ## {#HDF5_par_}
	 *
	 * Because the function **load** works independently from the number
	 * of processors the file created with save can be used also to post-process
	 * the saved data.
	 * Once loaded the particles on the distributed vector **vd2** we can
	 * use **vd2** to post-process the data. In this case we calculate
	 * the magnitude of the property 0 and we write it to a vtk file.
	 *
	 * \snippet Vector/1_HDF5_save_load/main.cpp hdf5_post_process
	 *
	 */

	//! \cond [hdf5_post_process] \endcond

	vector_dist<3,float, aggregate<float> > vd3(vd.getDecomposition(),0);

	auto it2 = vd2.getDomainIterator();

	while (it2.isNext())
	{
		auto p = it2.get();

		float magn_vort = sqrt(vd2.getProp<0>(p)[0]*vd2.getProp<0>(p)[0] +
				 	 	 	   vd2.getProp<0>(p)[1]*vd2.getProp<0>(p)[1] +
							   vd2.getProp<0>(p)[2]*vd2.getProp<0>(p)[2]);

		if (magn_vort < 0.1)
		{
			++it2;
			continue;
		}

		vd3.add();

		vd3.getLastPos()[0] = vd2.getPos(p)[0];
		vd3.getLastPos()[1] = vd2.getPos(p)[1];
		vd3.getLastPos()[2] = vd2.getPos(p)[2];

		vd3.template getLastProp<0>() = magn_vort;

		++it2;
	}

	// We write the vtk file out from vd3
	vd3.write("output", VTK_WRITER | FORMAT_BINARY );
	vd3.save("particles_post_process_2");

	//! \cond [hdf5_post_process] \endcond

	/*!
	 * \page Vector_1_HDF5 HDF5 save and load
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
