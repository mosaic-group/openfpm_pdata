//
// Created by jstark on 2020-05-18.
//
/**
 * @file example_sussman_images_2D/main.cpp
 * @page example_sussman_images_2D Images 2D
 *
 * # Load the geometrical object from a binary 2D image #
 * In this example, we will:
 * * Read a 2D binary image from a binary file (for a 3D image volume see @ref example_sussman_images_3D)
 * * Build a 2D cartesian OpenFPM grid with same dimensions as the image (1 particle for each pixel in x and y) or
 *   refined by arbitrary factor in dimension of choice (e.g. to get a isotropic grid)
 * * Assign pixel value to a grid node property
 * * Run the Sussman Redistancing and get the narrow band around the interface (see also: @ref example_sussman_circle
 *   and example_sussman_sphere)
 *
 * Output:
 * print on promt : Iteration, Change, Residual (see: #DistFromSol::change, #DistFromSol::residual)
 * writes vtk and hdf5 files of:
 * 1.) 2D grid with geometrical object pre-redistancing and post-redistancing (Phi_0 and Phi_SDF, respectively)
 * 2.) particles on narrow band around interface.
 *
 * ## Visualization of example output in Paraview ##
 * @htmlonly
 * <img src="http://openfpm.mpi-cbg.de/web/images/examples/sussman_redistancing/example_sussman_images_2D_paraview.png" width="1024px"/>
 * @endhtmlonly
 **/

/**
 * @page example_sussman_images_2D Images 2D
 *
 * ## Include ## {#e2d_img_include}
 *
 * These are the header files that we need to include:
 *
 * @snippet example/Numerics/Sussman_redistancing/example_sussman_images_2D/main.cpp Include
 *
 */
//! @cond [Include] @endcond

// Include redistancing header files
#include "util/PathsAndFiles.hpp"
#include "level_set/redistancing_Sussman/RedistancingSussman.hpp"
#include "level_set/redistancing_Sussman/NarrowBand.hpp"
#include "RawReader/InitGridWithPixel.hpp"
//! @cond [Include] @endcond


/**
 * @page example_sussman_images_2D Images 2D
 *
 * ## Initialization, indices and output folder ## {#e2d_img_init}
 *
 * Again we start with
 * * Initializing OpenFPM
 * * Setting the output path and creating an output folder
 * This time, we also set the input path and name of the binary image that we want to load onto the grid. For this
 * example we provide 3 simple example images. The original (e.g. tiff) image has been converted into -1 / +1 values.
 * A jupyter notebook that does this can be found here:
 * @htmlonly
 * <a href='http://openfpm.mpi-cbg.de/web/images/examples/sussman_redistancing/image2binary_dolphin.ipynb'>image2binary_dolphin.ipynb</a>.
 * @endhtmlonly
 * Optionally, we can define the grid dimensionality and some indices for better code readability later on.
 * * \p x: First dimension
 * * \p y: Second dimension
 * * \p Phi_0_grid: Index of property that stores the initial level-set-function
 * * \p Phi_SDF_grid: Index of property where the redistancing result should be written to
 *
 * @snippet example/Numerics/Sussman_redistancing/example_sussman_images_2D/main.cpp Initialization
 *
 */
//! @cond [Initialization] @endcond

int main(int argc, char* argv[])
{
//	std::string image_name = "circle";
//  std::string image_name = "triangle";
//  std::string image_name = "face";
	std::string image_name = "dolphin";


	//	initialize library
	openfpm_init(&argc, &argv);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Set current working directory, define output paths and create folders where output will be saved
	std::string cwd                     = get_cwd();
	const std::string path_output       = cwd + "/output_" + image_name;
//	const std::string path_output       = cwd + "/output_face";

	create_directory_if_not_exist(path_output);
	
	//////////////////////////////////////////////////////////////////////////////////////////////
	// Now we set the input paths. We need a binary file with the pixel values and a csv file with the
	// size of the stack (in #pixels / dimension)
	const std::string path_input        ="input/";
	const std::string path_to_image     = path_input + image_name + ".bin";
	const std::string path_to_size      = path_input + "size_" + image_name + ".csv";

	/*
	 * in case of 2D (single image):
	 */
	const unsigned int grid_dim         = 2;
	// some indices
	const size_t x                      = 0;
	const size_t y                      = 1;
	
	const size_t Phi_0_grid             = 0;
	const size_t Phi_SDF_grid           = 1;
	
	const size_t Phi_SDF_vd             = 0;
	const size_t Phi_grad_vd            = 1;
	const size_t Phi_magnOfGrad_vd      = 2;
	//! @cond [Initialization] @endcond
	
	
	/**
	 * @page example_sussman_images_2D Images 2D
	 *
	 * ## Set refinement factor ## {#e2d_img_refine}
	 *
	 * If we want that the grid has a different resolution as the image, we can set a refinement factor. In the
	 * refinement array, we define by which factor the grid resolution should be changed w.r.t. the image resolution.
	 * This can be useful for example when you want to have an isotropic grid but the underlying image is
	 * anisotropic (as it often happens for microcsopic z-stacks rather than 2D images).
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_images_2D/main.cpp Refinement
	 *
	 */
	//! @cond [Refinement] @endcond
	const double refinement []          = {0.5, 0.5};   // (2D) factor by which grid should be finer / coarser as
	// underlying image in each dimension (e.g. to get isotropic grid from anisotropic image resolution)
	//! @cond [Refinement] @endcond
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	read the stack size (number of pixel values per dimension) from a binary file
	// alternatively, you can directly define the stack size as: std::vector<size_t> stack_size {#pixels in x, #pixels in y}
	/**
	 * @page example_sussman_images_2D Images 2D
	 *
	 * ## Get size from image_size.csv ## {#e2d_img_size}
	 *
	 * Here we read the image size (number of pixel values per dimension) from a csv file. Alternatively, you can
	 * directly define the image size as: \p std::vector<size_t> \p stack_size \p {#pixels \p in \p x, \p #pixels \p in
	 * \p y}.
	 * We use this image size and the refinement factor to set the grid size \p sz.
	 *
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_images_2D/main.cpp Size
	 *
	 */
	//! @cond [Size] @endcond
	std::vector<size_t> stack_size = get_size(path_to_size);
	auto & v_cl = create_vcluster();
	if (v_cl.rank() == 0)
	{
		for(std::vector<int>::size_type i = 0; i != stack_size.size(); i++)
		{
			std::cout << "#Pixel in dimension " << i << " = " << stack_size[i] << std::endl;
			std::cout << "original refinement in dimension " << i << " = " << refinement[i] << std::endl;
		}
	}
	
	//  Array with grid dimensions after refinement. This size-array will be used for the grid creation.
	const size_t sz[grid_dim] = {(size_t)std::round(stack_size[x] * refinement[x]), (size_t)std::round(stack_size[y] * refinement[y])}; // 2D
	//! @cond [Size] @endcond
	// Here we create a 2D grid that stores 2 properties:
	// Prop1: store the initial Phi;
	// Prop2: here the re-initialized Phi (signed distance function) will be written to in the re-distancing step
	/**
	 * @page example_sussman_images_2D Images 2D
	 *
	 * ## Run Sussman redistancing and get narrow band ## {#e2d_img_redist}
	 *
	 * Once we have loaded the geometrical object from the 2D image onto the grid, we can perform Sussman
	 * redistancing and get the narrow band the same way as it is explained in detail here: @ref
	 * example_sussman_circle and here: @ref example_sussman_sphere.
	 *
	 *
	 * @snippet example/Numerics/Sussman_redistancing/example_sussman_images_2D/main.cpp Redistancing
	 *
	 */
	//! @cond [Redistancing] @endcond
	Box<grid_dim, double> box({0.0, 0.0}, {5.0, 5.0}); // 2D
	
	Ghost<grid_dim, long int> ghost(0);
	typedef aggregate<double, double> props;
	typedef grid_dist_id<grid_dim, double, props > grid_in_type;
	grid_in_type g_dist(sz, box, ghost);
	g_dist.setPropNames({"Phi_0", "Phi_SDF"});
	
	// Now we can initialize the grid with the pixel values from the image stack
	load_pixel_onto_grid<Phi_0_grid>(g_dist, path_to_image, stack_size);
	g_dist.write(path_output + "/grid_from_images_initial", FORMAT_BINARY); // Save the initial grid as vtk file


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Now we want to convert the initial Phi into a signed distance function (SDF) with magnitude of gradient = 1
	// For the initial re-distancing we use the Sussman method
	// 1.) Set some redistancing options
	Redist_options redist_options;
	redist_options.min_iter                             = 100;      // min. number of iterations before steady state in narrow band will be checked (default: 100)
	redist_options.max_iter                             = 10000;    // max. number of iterations you want to run the
																	// redistancing, even if steady state might not yet have been reached (default: 1e6)
	
	redist_options.convTolChange.value                  = 1e-6;     // convolution tolerance for the normalized total change of Phi in the narrow band between two consecutive iteratins (default: 1e-6)
	redist_options.convTolChange.check                  = true;     // define here which of the convergence criteria above should be used. If both are true, termination only occurs when both are fulfilled or when iter > max_iter
	redist_options.convTolResidual.value                = 1e-1;     // convolution tolerance for the normalized total residual of Phi versus the SDF (aka the error). Don't choose too small to not run forever. (default: 1e-1)
	redist_options.convTolResidual.check                = false;    // (default: false)
	
	redist_options.interval_check_convergence           = 1;        // interval of #iterations at which convergence is checked (default: 100)
	redist_options.width_NB_in_grid_points              = 6;        // width of narrow band in number of grid points.
	// Must be at least 4, in order to have at least 2 grid points on each side of the interface. (default: 4)
	redist_options.print_current_iterChangeResidual     = true;     // if true, prints out every current iteration + corresponding change from the previous iteration + residual from SDF (default: false)
	redist_options.print_steadyState_iter               = true;     // if true, prints out the final iteration number when steady state was reached + final change + residual (default: true)
	
	RedistancingSussman<grid_in_type> redist_obj(g_dist, redist_options);   // Instantiation of Sussman-redistancing class
//	std::cout << "dt = " << redist_obj.get_time_step() << std::endl;
	// Run the redistancing. in the <> brackets provide property-index where 1.) your initial Phi is stored and 2.) where the resulting SDF should be written to.
	redist_obj.run_redistancing<Phi_0_grid, Phi_SDF_grid>();
	
	g_dist.write(path_output + "/grid_images_post_redistancing", FORMAT_BINARY); // Save the result as vtk file
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Get narrow band: Place particles on interface (narrow band width e.g. 2 grid points on each side of the interface)
	size_t bc[grid_dim] = {PERIODIC, PERIODIC};
	// Create an empty vector to which narrow-band particles will be added. You can choose, how many properties you want.
	// Minimum is 1 property, to which the Phi_SDF can be written
	// In this example we chose 3 properties. The 1st for the Phi_SDF, the 2nd for the gradient of phi and the 3rd for
	// the magnitude of the gradient
	typedef aggregate<double, double[grid_dim], double> props_nb;
	typedef vector_dist<grid_dim, double, props_nb> vd_type;
        Ghost<grid_dim, double> ghost_vd(0);
        vd_type vd_narrow_band(0, box, bc, ghost_vd);
	vd_narrow_band.setPropNames({"Phi_SDF", "Phi_grad", "Phi_magnOfGrad"});
	
	NarrowBand<grid_in_type> narrowBand(g_dist, redist_options.width_NB_in_grid_points); // Instantiation of NarrowBand class
	
	// Get the narrow band. You can decide, if you only want the Phi_SDF saved to your particles or
	// if you also want the gradients or gradients and magnitude of gradient.
	// The function will know depending on how many property-indices you provide in the <> brackets.
	// First property-index must always be the index of the SDF on the grid!
	// E.g.: The following line would only write only the Phi_SDF from the grid to your narrow-band particles
	// narrowBand.get_narrow_band<Phi_SDF_grid, Phi_SDF_vd>(g_dist, vd_narrow_band);
	// Whereas this will give you the gradients and magnOfGrad of Phi as well:
	narrowBand.get_narrow_band<Phi_SDF_grid, Phi_SDF_vd, Phi_grad_vd, Phi_magnOfGrad_vd>(g_dist, vd_narrow_band);
	
	vd_narrow_band.write(path_output + "/vd_narrow_band_images", FORMAT_BINARY); // Save particles as vtk file

	openfpm_finalize(); // Finalize openFPM library
	return 0;
}
//! @cond [Redistancing] @endcond

/**
 * @page example_sussman_images_2D Images 2D
 *
 * ## Full code ## {#e2d_img_full}
 *
 * @include example/Numerics/Sussman_redistancing/example_sussman_images_2D/main.cpp
 */









