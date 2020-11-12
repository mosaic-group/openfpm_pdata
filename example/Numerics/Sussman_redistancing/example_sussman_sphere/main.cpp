//
// Created by jstark on 2020-05-18.
//
/**
 * @file example_sussman_sphere/main.cpp
 * @page example_sussman_sphere Sphere 3D
 *
 * # Get the signed distance function for a 3D sphere via Sussman Redistancing #
 *
 * In this example, we perform Sussman redistancing on a filled 3D sphere. Here, the sphere is constructed via
 * equation. For image based reconstruction and redistancing see @ref example_sussman_images_3D.
 *
 * A filled sphere with center a, b, c can be constructed by the following equation:
 *
 *
 * @f[ (x-a)^2 + (y-b)^2 + (z-c)^2 <= r^2 @f]
 *
 * We will create a filled sphere on a 3D grid, where the sphere is represented by a -1/+1 step function
 * (indicator function) as:
 *
 * @f[ \phi_{\text{indicator}} = \begin{cases}
 *    +1 & \text{point lies inside the object} \\
 *     0 & \text{object surface} \\
 *    -1 & \text{point lies outside the object} \\
 * \end{cases} @f]
 *
 * By Sussman redistancing, this indicator function is transformed into a signed distance function:
 *
 * @f[ \phi_{\text{SDF}} = \begin{cases}
 *    +d & \text{orthogonal distance to the closest surface of a point lying inside the object} \\
 *     0 & \text{object surface} \\
 *    -d & \text{orthogonal distance to the closest surface of a point lying outside the object} \\
 * \end{cases} @f]
 *
 * Once we have received the Phi_SDF from the redistancing, particles can be placed on narrow band around the interface.
 *
 * * Creates filled 3D sphere with -1/+1 indicator function
 * * Runs Sussman redistancing (see @ref RedistancingSussman.hpp)
 * * Places particles on narrow band around interface
 *
 * Output:
 * print on promt : Iteration, Change, Residual (see: #DistFromSol::change, #DistFromSol::residual)
 * writes vtk and hdf5 files of:
 * 1.) 3D grid with sphere pre-redistancing and post-redistancing (Phi_0 and Phi_SDF, respectively)
 * 2.) particles on narrow band around interface
 *
 * ## Visualization of example output in Paraview ##
 * @htmlonly
 * <img src="http://openfpm.mpi-cbg.de/web/images/examples/sussman_redistancing/example_sussman_sphere_paraview.png" width="1024px"/>
 * @endhtmlonly
 *
 */

/**
 * @page example_sussman_sphere Sphere 3D
 *
 * ## Include ## {#e3d_s_include}
 *
 * These are the header files that we need to include:
 *
 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Include
 *
 */
//! @cond [Include] @endcond

// Include level-set-method related header files
#include "util/PathsAndFiles.hpp"
#include "level_set/redistancing_Sussman/RedistancingSussman.hpp"
#include "level_set/redistancing_Sussman/NarrowBand.hpp"
#include "Draw/DrawSphere.hpp"
//! @cond [Include] @endcond

/**
 * @page example_sussman_sphere Sphere 3D
 *
 * ## Initialization and output folder ## {#e2d_s_init}
 *
 * We start with
 * * Initializing OpenFPM
 * * Setting the output path and creating an output folder
 *
 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Initialization and output folder
 *
 */
//! @cond [Initialization and output folder] @endcond
int main(int argc, char* argv[])
{
	//	initialize library
	openfpm_init(&argc, &argv);
	// Set current working directory, define output paths and create folders where output will be saved
	std::string cwd                     = get_cwd();
	const std::string path_output       = cwd + "/output_sphere";
	create_directory_if_not_exist(path_output);
	//! @cond [Initialization and output folder] @endcond
	
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Indices for the grid ## {#e3d_s_indices}
	 *
	 * Optionally, we can define the grid dimensionality and some indices for better code readability later on.
	 * * \p x: First dimension
	 * * \p y: Second dimension
	 * * \p z: Third dimension
	 * * \p Phi_0_grid: Index of property that stores the initial level-set-function
	 * * \p Phi_SDF_grid: Index of property where the redistancing result should be written to
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Indices grid
	 *
	 */
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//! @cond [Indices grid] @endcond	
	// grid dimension
	constexpr size_t grid_dim = 3;
	// Some indices for the grid / grid properties
	const size_t x                      = 0;
	const size_t y                      = 1;
	const size_t z                      = 2;
	
	const size_t Phi_0_grid             = 0;
	const size_t Phi_SDF_grid           = 1;
	//! @cond [Indices grid] @endcond
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Create the grid ## {#e3d_s_grid}
	 *
	 * In order to get a grid, we
	 * * Define the grid size in terms of number of grid points per dimension
	 * * Create a Box that defines our domain
	 * * Create a Ghost object that will define the extension of the ghost part
	 * * Create a 3D grid with two properties of type double, one for the pre-redistancing Phi_initial and one, where
	 *   the post-redistancing Phi_SDF should be written to.
	 * * Set some property names (optionally. These names show up when opening the grid vtk in Paraview.)
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Grid creation
	 *
	 */
	//! @cond [Grid creation] @endcond
	// Here we create a 3D grid that stores two properties:
	// Prop1: store the initial Phi;
	// Prop2: here the re-initialized Phi (signed distance function) will be written to in the re-distancing step
	const size_t sz[grid_dim] = {128, 128, 128};
	Box<grid_dim, double> box({0.0, 0.0, 0.0}, {10.0, 10.0, 10.0});
	Ghost<grid_dim, long int> ghost(0);
	typedef aggregate<double, double> props;
	typedef grid_dist_id<grid_dim, double, props > grid_in_type;
	grid_in_type g_dist(sz, box, ghost);
	g_dist.setPropNames({"Phi_0", "Phi_SDF"});
	//! @cond [Grid creation] @endcond
	
	
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Get filled sphere on the grid ## {#e3d_s_getsphere}
	 *
	 *
	 * On the grid that we have just created, we can now initialize Phi_0 as a filled sphere of defined radius.
	 * The center of the sphere is passed as x_center and y_center (see @ref Sphere.hpp). Phi_0 will then be:
	 *
	 * @f[ \phi_{\text{indicator}} = \begin{cases}
     *    +1 & \text{point lies inside the object} \\
     *     0 & \text{object surface} \\
	 *    -1 & \text{point lies outside the object} \\
     * \end{cases} @f]
	 *
	 * Optionally, we can save this initial grid as a vtk file, if we want to visualize and check it in Paraview.
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Get sphere
	 */
	//! @cond [Get sphere] @endcond
	// Now we initialize the grid with a filled sphere. Outside the sphere, the value of Phi_0 will be -1, inside +1.
	double radius = 1.0; // Radius of the sphere
	init_grid_with_sphere<Phi_0_grid>(g_dist, radius, 5, 5, 5); // Initialize sphere onto grid, centered at (5, 5, 5)

	g_dist.write(path_output + "/grid_initial_sphere_preRedistancing_radius" + std::to_string((int)radius) , FORMAT_BINARY); // Save the sphere as vtk file
	//! @cond [Get sphere] @endcond



    /**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Set the redistancing options ## {#e3d_s_redistoptions}
	 *
	 * For the redistancing, we can choose some options. These options will then be passed bundled as a structure to
	 * the redistancing function. Setting these options is optional, since they all have a Default value as well. In
	 * particular the following options can be set by the user:
	 * * \p min_iter: Minimum number of iterations before steady state in narrow band will be checked (Default: 100).
	 * * \p max_iter: Maximum number of iterations you want to run the redistancing, even if steady state might not yet
	 *                have been reached (Default: 1e6).
	 * * \p convTolChange.value: Convolution tolerance for the normalized total change of Phi in the narrow band between
	 *                           two consecutive iterations (Default: 1e-6).
	 * * \p convTolChange.check: Set true, if you want to use the normalized total change between two iterations as
	 *                           measure of how close you are to the steady state solution. Redistancing will then stop
	 *                           if convTolChange.value is reached or if the current iteration is bigger than max_iter.
	 * * \p convTolResidual.value: Convolution tolerance for the residual, that is abs(magnitude gradient of phi - 1) of
	 *                             Phi in the narrow band (Default 1e-1).
	 * * \p convTolResidual.check: Set true, if you want to use the residual of the current iteration as measure of how
	 *                             close you are to the steady state solution. Redistancing will then stop if
	 *                             convTolResidual.value is reached or if the current iteration is bigger than max_iter.
	 * * \p interval_check_convergence: Interval of #iterations at which convergence to steady state is checked
	 *                                  (Default: 100).
	 * * \p width_NB_in_grid_points: Width of narrow band in number of grid points. Must be at least 4, in order to 
	 *                               have at least 2 grid points on each side of the interface. Is automatically set 
	 *                               to 4, if a value smaller than 4 is chosen (Default: 4).
	 * * \p print_current_iterChangeResidual: If true, the number of the current iteration, the corresponding change
	 *                                        w.r.t the previous iteration and the residual is printed (Default: false).
	 * * \p print_steadyState_iter: If true, the number of the steady-state-iteration, the corresponding change
	 *                              w.r.t the previous iteration and the residual is printed (Default: false).
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Redistancing options
	 */
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//! @cond [Redistancing options] @endcond
	// Now we want to convert the initial Phi into a signed distance function (SDF) with magnitude of gradient = 1.
	// For the initial re-distancing we use the Sussman method. First of all, we can set some redistancing options.
	Redist_options redist_options;
	redist_options.min_iter                             = 100;
	redist_options.max_iter                             = 1000;
	
	redist_options.convTolChange.value                  = 1e-12;
	redist_options.convTolChange.check                  = true;
	redist_options.convTolResidual.value                = 1e-1;
	redist_options.convTolResidual.check                = false;
	
	redist_options.interval_check_convergence           = 1;
	redist_options.width_NB_in_grid_points              = 6;
	redist_options.print_current_iterChangeResidual     = true;
	redist_options.print_steadyState_iter               = true;
	//! @cond [Redistancing options] @endcond
	
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Run the redistancing ## {#e3d_s_runredist}
	 *
	 * Now, we can instantiate the Sussman-redistancing class using the redist_options we have just set (or use the
	 * defaults) and run the redistancing. If we are interested in the time-step that was automatically computed
	 * internally according to the CFL-condition, we can access it with #RedistancingSussman::get_time_step(). If, for
	 * whatever reason, we wish to change this timestep, we can do that with #RedistancingSussman::set_user_time_step
	 * (T dt). However, be careful when setting your own timestep: if chosen too large, the CFL-condition may be no
	 * longer fulfilled and the method can become unstable. For the #run_redistancing we need to provide two template
	 * parameters. These are the indices of the respective property that: 1.) contains the initial Phi, 2.) should
	 * contain the output of the redistancing. You can use the same index twice, if you want that your input will be
	 * overwritten by the output. This time, it makes sense, to save the output grid, since this is already our grid
	 * containing the signed distance function in Prop. 2. The vtk-file can be opened in Paraview. If we want, we can
	 * further save the result as hdf5 file that can be reloaded onto an openFPM grid.
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Run redistancing
	 */
	//! @cond [Run redistancing] @endcond
	RedistancingSussman<grid_in_type> redist_obj(g_dist, redist_options);   // Instantiation of Sussman-redistancing class
    //	std::cout << "dt = " << redist_obj.get_time_step() << std::endl;
    // Run the redistancing. in the <> brackets provide property-index where 1.) your initial Phi is stored and 2.) 
    // where the resulting SDF should be written to.
	redist_obj.run_redistancing<Phi_0_grid, Phi_SDF_grid>();
	
	g_dist.write(path_output + "/grid_sphere_postRedistancing", FORMAT_BINARY); // Save the result as vtk file
	g_dist.save(path_output + "/grid_sphere_postRedistancing" + ".bin"); // Save the result as hdf5 file that can be 
	// reloaded onto an openFPM grid
	//! @cond [Run redistancing] @endcond
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Initialize empty vector for narrow band particles ## {#e3d_s_initnbvector}
	 *
	 * Next, we want to place particles in a narrow band around the interface. For this, we need to create an empty
	 * vector in the same domain as our grid from before. We therefore re-use the \p box and \p ghost that we defined
	 * during grid creation. This time, we have to set the boundary conditions (see general OpenFPM vector_dist
	 * documentation for details.) We also have to choose now, how many properties our particle should store. Minimum is
	 * one, but you can have particles with arbitrary many properties, depending on what you want to use them for
	 * later on. Here, we exemplarily define 3 properties.
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Initialize narrow band
	 */
	//! @cond [Initialize narrow band] @endcond
	//	Get narrow band: Place particles on interface (narrow band width e.g. 2 grid points on each side of the interface)
	size_t bc[grid_dim] = {PERIODIC, PERIODIC, PERIODIC};
	// Create an empty vector to which narrow-band particles will be added. You can choose, how many properties you want.
	// Minimum is 1 property, to which the Phi_SDF can be written
	// In this example we chose 3 properties. The 1st for the Phi_SDF, the 2nd for the gradient of phi and the 3rd for
	// the magnitude of the gradient
	typedef aggregate<double, double[grid_dim], double> props_nb;
	typedef vector_dist<grid_dim, double, props_nb> vd_type;
        Ghost<grid_dim, double> ghost_vd(0);
        vd_type vd_narrow_band(0, box, bc, ghost_vd);
	vd_narrow_band.setPropNames({"Phi_SDF", "Phi_grad", "Phi_magnOfGrad"});
	//! @cond [Initialize narrow band] @endcond
	
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Instantiate narrow band ## {#e3d_s_instnarrowband}
	 *
	 * Before we fill the empty vector that we just created with particles, we define how thick the narrow band
	 * around the interface should be.
	 * This is independent of the redist_options.width_NB_in_grid_points that we chose before. However, it makes
	 * sense to use a width during redistancing which is at least as large as the width that we want to place the
	 * particles in. The reason is that convergence criteria were checked for that part of the grid only, which
	 * belonged to the narrow band.
	 *
	 * Depending on the type of the variable which you define, the width can be either set in terms of number of
	 * grid points (size_t), physical width (double) or extension of narrow band as physical width inside of object
	 * and outside the object (double, double).
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Instantiate narrow band
	 */
	//! @cond [Instantiate narrow band] @endcond
	size_t thickness_of_narrowBand_in_grid_points = 6;
	NarrowBand<grid_in_type> narrowBand(g_dist, redist_options.width_NB_in_grid_points); // Instantiation of NarrowBand class

//! @cond [Instantiate narrow band] @endcond
	
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Indices for the narrow band vector ## {#e3d_s_indices}
	 *
	 * Again, we can define some indices for better code readability. This is just an example, you may want to choose
	 * different names and have a different number of properties thus different number of indices.
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Indices narrow band
	 *
	 */
	//////////////////////////////////////////////////////////////////////////////////////////////
	//! @cond [Indices narrow band] @endcond
	// Some indices
	const size_t Phi_SDF_vd             = 0;
	const size_t Phi_grad_vd            = 1;
	const size_t Phi_magnOfGrad_vd      = 2;
	//! @cond [Indices narrow band] @endcond
	
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Fill the vector with particles placed within a narrow band around the interface ## {#e3d_s_getnarrowband}
	 *
	 * Now, we finally fill the vector with particles placed in a narrow band around the interface. The first
	 * template parameter that we pass, must always be the index of the grid property that contains the Phi_SDF.
	 *
	 * Then you can decide between the following options, what you want to copy from the grid to the particles:
	 * * Only copy the Phi_SDF: \p narrowBand.get_narrow_band<Phi_SDF_grid, Phi_SDF_vd>(g_dist, vd_narrow_band);
	 * * Copy Phi_SDF as well as the first order upwind gradient of Phi_SDF: \p narrowBand
	 *   .get_narrow_band<Phi_SDF_grid, Phi_SDF_vd, Phi_grad_vd>(g_dist, vd_narrow_band);
	 * * Copy Phi_SDF as well as the first order upwind gradient of Phi_SDF and the corresponding gradient magnitude:
	 *   \p narrowBand.get_narrow_band<Phi_SDF_grid, Phi_SDF_vd, Phi_grad_vd, Phi_magnOfGrad_vd>(g_dist,
	 *   vd_narrow_band);
	 * * Copy an arbitrary property: \p narrowBand.get_narrow_band_copy_specific_property<Phi_SDF_grid,
	 *   Prop_index_grid, Prop_index_vector>(g_dist, vd_narrow_band);
	 *
	 * We save the particles in a vtk file (open in Paraview) and as hdf5 file (can be loaded back on particles).
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Get narrow band
	 */
	// Get the narrow band. You can decide, if you only want the Phi_SDF saved to your particles or
	// if you also want the gradients or gradients and magnitude of gradient.
	// The function will know depending on how many property-indices you provide in the <> brackets.
	// First property-index must always be the index of the SDF on the grid!
	// E.g.: The following line would only write only the Phi_SDF from the grid to your narrow-band particles
	// narrowBand.get_narrow_band<Phi_SDF_grid, Phi_SDF_vd>(g_dist, vd_narrow_band);
	// Whereas this will give you the gradients and magnOfGrad of Phi as well:
	//! @cond [Get narrow band] @endcond
	narrowBand.get_narrow_band<Phi_SDF_grid, Phi_SDF_vd, Phi_grad_vd, Phi_magnOfGrad_vd>(g_dist, vd_narrow_band);
	
	vd_narrow_band.write(path_output + "/vd_narrow_band_sphere", FORMAT_BINARY); // Save particles as vtk file
	vd_narrow_band.save(path_output + "/vd_narrow_band_sphere.bin"); // Save particles as hdf5 file -> can be reloaded as particles
	//! @cond [Get narrow band] @endcond
	
	/**
	 * @page example_sussman_sphere Sphere 3D
	 *
	 * ## Terminate ## {#e3d_s_terminate}
	 *
	 * We end with terminating OpenFPM
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp Terminate
	 */
	//! @cond [Terminate] @endcond
	openfpm_finalize(); // Finalize openFPM library
	return 0;
}
//! @cond [Terminate] @endcond

/**
 * @page example_sussman_sphere Sphere 3D
 *
 * ## Full code ## {#e3d_s_full}
 *
 * @include example/Numerics/Sussman_redistancing/example_sussman_sphere/main.cpp
 */




































