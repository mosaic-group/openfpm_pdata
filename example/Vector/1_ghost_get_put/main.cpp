/*!
 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
 *
 *
 * [TOC]
 *
 *
 * # Ghost get and ghost put # {#ghost_get_put}
 *
 *
 * This example shows more in details the functionalities of **ghost_get** and **ghost_put** for a distributed vector.
 *
 * ## Inclusion ## {#e1_v_inclusion}
 *
 * We activate the vector_dist functionalities
 *
 * \snippet Vector/1_ghost_get_put/main.cpp inclusion
 *
 * \see \ref e0_v_inclusion
 *
 */

//! \cond [inclusion] \endcond
#include "Vector/vector_dist.hpp"
//! \cond [inclusion] \endcond

int main(int argc, char* argv[])
{

	/*!
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## Initialization ## {#e1_s_init}
	 *
	 *  Here we
	 *  * Initialize the library
	 *  * we create a Box that define our domain
	 *  * An array that define our boundary conditions
	 *  * A Ghost object that will define the extension of the ghost part in physical units
	 *
	 * \see \ref e0_s_init
	 *
	 * \snippet Vector/1_ghost_get_put/main.cpp Initialization and parameters
	 *
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
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## Vector instantiation ## {#e1_v_init}
	 *
	 * Here we are creating a distributed vector following the previous example
	 *
	 * \see \ref e0_s_vector_inst
	 *
	 * \snippet Vector/1_ghost_get_put/main.cpp vector instantiation
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
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## Assign values ## {#e1_s_assign_pos}
	 *
	 * We get an iterator and we iterate across the 4096 particles where we define their positions and properties
	 *
	 * \see \ref e0_s_assign_pos
	 *
	 *  \snippet Vector/1_ghost_get_put/main.cpp assign position
	 *
	 */

	//! \cond [assign position] \endcond

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		// Particle p
		auto p = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(p)[0] = (float)rand() / RAND_MAX;

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(p)[1] = (float)rand() / RAND_MAX;

        // the scalar property
		vd.template getProp<scalar>(p) = vd.getPos(p)[0];

		// The vector
		vd.template getProp<vector>(p)[0] = 1.0;
		vd.template getProp<vector>(p)[1] = 2.0;
		vd.template getProp<vector>(p)[2] = 3.0;

		// The tensor
		vd.template getProp<tensor>(p)[0][0] = 1.0;
		vd.template getProp<tensor>(p)[0][1] = 2.0;
		vd.template getProp<tensor>(p)[0][2] = 3.0;
		vd.template getProp<tensor>(p)[1][0] = 4.0;
		vd.template getProp<tensor>(p)[1][1] = 5.0;
		vd.template getProp<tensor>(p)[1][2] = 6.0;
		vd.template getProp<tensor>(p)[2][0] = 7.0;
		vd.template getProp<tensor>(p)[2][1] = 8.0;
		vd.template getProp<tensor>(p)[2][2] = 9.0;

		// next particle
		++it;
	}

	//! \cond [assign position] \endcond

	/*!
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## Mapping particles ## {#e1_s_map}
	 *
	 * We redistribute the particles according to the underline decomposition
	 *
	 * \snippet Vector/1_ghost_get_put/main.cpp map
	 *
	 *
	 */

	//! \cond [map] \endcond

	vd.map();

	//! \cond [map] \endcond

	/*!
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## ghost_get ## {#e1_v_assign_prop}
	 *
	 * Here we synchronize the ghosts in the standard way, updating the scalar and
	 * the tensor property.
	 *
	 * \warning Every **ghost_get** by default reset the status of the ghost so the previous information is lost
	 *
	 * In the 2 images below we see what happen when we do a standard **ghost_get**
	 * Before and after. The blue arrows in the first image indicate the vector field
	 * for the real particles. In the second image instead the red arrow indicate the
	 * vector field for the real particle. The blue arrow indicate the ghosts. We can
	 * note that the blue arrow does not contain the correct vector. The reason is that
	 * when we used **ghost_get** we synchronized the scalar, and the tensor, but not the vector.
	 *
	 * \see \ref e1_part_ghost
	 *
	 * \snippet Vector/1_ghost_get_put/main.cpp ghost_get_some_prop
	 *
	 * \htmlonly
	 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/before_ghost_get.jpg"/>
	 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/after_ghost_get.jpg"/>
	 * \endhtmlonly
	 *
	 * ## So ... how I have to put these ghost_get ##
	 *
	 * The first thing to do is to place the ghost in a way that the program work
	 * in parallel for sure. In order to do this we can do the following reasoning:
	 * If we have a loop over particles we distinguish two type of loops:
	 *
	 *  * A loop that iterate over particles
	 *  * A loop that iterate over particles and neighborhood particles
	 *
	 *
	 * If the loop is of the first type (you do not loop over the neighborhood particles)
	 *  ghost_get is not necessary. If I am in the second case I need a ghost_get. The
	 *  second point is which property I have to synchronize ghost_get<...>(), or more
	 *  practically what I have to put in the ... . To answer this we have to check all
	 *  the properties that we use from the neighborhood particles and pass it to ghost_get
	 *  as a list. To summarize:
       \code{.unparsed}

        I am doing a simple loop over particles (1), or I am looping also over neighborhood particles (2)?
        For the case (1) the answer is "I do not need ghost_get". For the case (2) the answer is "I need ghost_get"

        if I am on the case (2) the second question is which parameters should I use ?
                  The answer is look at all vd.getProp<...>(b) where b is a neighborhood particle. All ... properties should appear in
                  ghost_get<...>()

       \endcode
     * This reasoning is always enough to have ghost_get function always placed correctly. For
     * more fine tuning look at the options below
	 *
	 */

	//! \cond [ghost_get_some_prop] \endcond

	// write in vtk to see the result before
	vd.write("before_ghost_get");

	// Here we synchronize the ghost get only the scalar and the tensor property
	vd.ghost_get<scalar,tensor>();

	// write in vtk to see the result after
	vd.write("after_ghost_get");

	//! \cond [ghost_get_some_prop] \endcond

	/*!
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## ghost_get KEEP_PROPERTIES ## {#e1_gg_options}
	 *
	 * As described before every ghost_get reset the status of the ghost. In particular
	 * all the informations are lost. So for example doing a **ghost_get<scalar>** followed
	 * by a **ghost_get<vector>** will not preserve the information of the **ghost_get<scalar>**.
	 *  This because in general **ghost_get**
	 * recompute the set of particles to send to the other processor based on their
	 * geometrical position. If the particles move between the **ghost_get<scalar>** and the
	 * the **ghost_get<vector>** the set of the particles sent can change between the 2 **ghost_get**.
	 * Merge the information between the two different set of particles would end up to be non trivial and would
	 * require to send additional information. This will make the communication inconveniently
	 * heavier. In the case the particles does not move on the other hand it would be
	 * possible to do a trivial merge of the information. The library is not able to
	 * detect automatically such cases. The information must be given explicitly.
	 * **KEEP_PROPERTIES** is an option that can be given to ghost_get to force to do
	 * a trivial merge in case the particles do not move. More in general such property is
	 * safe to be used in case that the particles move between **ghost_get**. What will
	 * happen is that the computation of the set of the particles to send will be fully skipped.
	 * OpenFPM will send the same set of particles from the previous ghost_get. Such functionality
	 *  is important in case of usage of Verlet-List with radius+skin. This functionality is
	 *   advanced and will be explained in a next example.
	 *
	 * \see \ref Vector_5_md_vl_sym
	 *
	 *  Because the option **KEEP_PROPERTIES** has the functionality to keep the properties
	 *  and skip the labelling. The option **KEEP_PROPERTIES** is also names **SKIP_LABELLING**.
	 * The two options are exactly equivalents.
	 *
	 * In the picture below we can see what happen when from the previous situation,
	 * we do a **ghost_get** with **KEEP_PROPERTIES** for the vector field. The vector change
	 * pointing to the right direction. On the other hand also the  scalar (The dot color)is kept and the information
	 * is not destroyed.
	 *
	 * \snippet Vector/1_ghost_get_put/main.cpp keep prop
	 *
	 * \htmlonly
	 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/after_ghost_get_kp.jpg"/>
	 * \endhtmlonly
	 *
	 *
	 */

	//! \cond [keep prop] \endcond

	vd.ghost_get<vector>(KEEP_PROPERTIES);
	
	// write in vtk to see the result before
	vd.write("after_ghost_get_kp");

	//! \cond [keep prop] \endcond

	/*!
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## ghost_get NO_POSITION ##
	 *
	 * **ghost_get** in general send automatically the information about position.
	 * If we are sure that the particles position did not change its position we
	 * can use the option NO_POSITION, to avoid to send the positional informations.
	 * With the only purpose to show what happens we shift the particle position of the
	 * ghost parts by one. We also force the vector to point to the opposite direction.
	 *
	 * \htmlonly
	 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/before_ghost_get_kp_no_pos.jpg"/>
	 * \endhtmlonly
	 *
	 * Doing a **ghost_get** with KEEP_PROPERTY and NO_POSITION, it will updates the information
	 * of the vectors, but the position will remain the same.
	 *
	 * \htmlonly
	 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/after_ghost_get_kp_no_pos.jpg"/>
	 * \endhtmlonly
	 *
	 * We can restore the positioninformation doing a **ghost_get** without NO_POSITION
	 *
	 * \htmlonly
	 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/after_ghost_get_kp_restore.jpg"/>
	 * \endhtmlonly
	 *
	 * \snippet Vector/1_ghost_get_put/main.cpp ghost_get no pos
	 *
	 *
	 *
	 */

	//! \cond [ghost_get no pos] \endcond

	it = vd.getGhostIterator();

	while (it.isNext())
	{
		// Particle p
		auto p = it.get();

		// we shift down he particles
		vd.getPos(p)[0] -= 1.0;

		// we shift
		vd.getPos(p)[1] -= 1.0;

		// The vector
		vd.template getProp<vector>(p)[0] = -1.0;
		vd.template getProp<vector>(p)[1] = -2.0;
		vd.template getProp<vector>(p)[2] = -3.0;

		// next particle
		++it;
	}

	vd.ghost_get<vector>(KEEP_PROPERTIES | NO_POSITION);

	// write in vtk to see the result before
	vd.write("after_ghost_get_kp_no_pos");

	vd.ghost_get<vector>(KEEP_PROPERTIES);

	vd.write("after_ghost_get_kp_no_pos_restore");

	//! \cond [ghost_get no pos] \endcond

	/*!
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## ghost_put ## {#ghost_put_e1}
	 *
	 * **ghost_put** is another particular function. It work similarly to **ghost_get**
	 *  but in the inverted direction. In particular it take the information from
	 *  the ghost particles and it send the information back to the real particles.
	 *  How to merge the information back to the real particles is defined by an operation.
	 *
	 * In this particular case we scatter back the information present in the ghost.
	 *
	 * \note **ghost_put** does not re-compute the set of
	 *       the particles to send. This mean that the set of particles communicated back from
	 *       ghost_put is given by the set of particles received by the last ghost_get. The
	 *       reason of this are similar to the **ghost_get**. If the set of particles is recomputed
	 *       there is no trivial way for the other processor to merge the information. This mean
	 *       that additional information must be sent and this is most of the time inconvenient.
	 *
	 * Once the particles are received back we add_ such contributions to the real particles.
	 * A typical real application of such functionality is in the case of symmetric interactions
	 *
	 * As we can see from the image below some of red particles near the ghost has a bigger arrow
	 *
	 * \htmlonly
	 * <img src="http://ppmcore.mpi-cbg.de/web/images/examples/after_ghost_put.jpg"/>
	 * \endhtmlonly
	 *
	 * \see \ref Vector_5_md_vl_sym
	 *
	 * \snippet Vector/1_ghost_get_put/main.cpp ghost_put
	 *
	 */

	//! \cond [ghost_put] \endcond

	vd.ghost_put<add_,vector>();

	// write in vtk to see the result before
	vd.write("after_ghost_put");

	//! \cond [ghost_put] \endcond

	/*!
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## Finalize ## {#finalize_e1_sim}
	 *
	 *  At the very end of the program we have always de-initialize the library
	 *
	 * \snippet Vector/1_ghost_get_put/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_1_ghost_get Vector 1 ghost get and ghost put
	 *
	 * ## Full code ## {#code_e1_sim}
	 *
	 * \include Vector/1_ghost_get_put/main.cpp
	 *
	 */
}
