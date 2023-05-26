//
// Created by jstark on 28.09.21.
//

#ifndef SPARSEGRID_DIFFUSIONSPACE_SPARSEGRID_HPP
#define SPARSEGRID_DIFFUSIONSPACE_SPARSEGRID_HPP

#include "math.h"

template <typename T>
static bool is_diffusionSpace(const T & phi_sdf, const T & lower_bound, const T & upper_bound)
{
	const T EPSILON = std::numeric_limits<T>::epsilon();
	const T _b_low = lower_bound + EPSILON;
	const T _b_up  = upper_bound  - EPSILON;
	return (phi_sdf > _b_low && phi_sdf < _b_up);
}

template <size_t PHI_SDF_ECS_FULL, size_t PHI_SDF_ECS_SPARSE, typename grid_type, typename sparse_in_type, typename T>
void get_diffusion_domain_sparse_grid(grid_type & grid,
									  sparse_in_type & sparse_grid,
									  const T b_low,
									  const T b_up)
{
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		if(is_diffusionSpace(grid.template get<PHI_SDF_ECS_FULL>(key), b_low, b_up))
		{
			sparse_grid.template insertFlush<PHI_SDF_ECS_SPARSE>(key) = grid.template get<PHI_SDF_ECS_FULL>(key);
		}
		++dom;
	}
	// Mapping?
	// Load Balancing?
}

#if 0
template <size_t PHI_SDF_ECS_FULL, size_t PHI_SDF_ECS_SPARSE, typename grid_type, typename grid_type_upscale, typename sparse_in_type, typename T>
void get_diffusion_domain_sparse_grid_upscale(grid_type & grid,
											  grid_type_upscale & grid_upscale,
									          sparse_in_type & sparse_grid,
									          const int upscale_factor, // total upscale factor that is product all dimensions, e.g. 2^3 = 8
									          const T b_low,
									          const T b_up)
{
	auto dom = grid.getDomainIterator();
	auto dom_upscale = grid_upscale.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		auto key_upscale = dom_upscale.get();

		for(int i = 0, i < upscale_factor, ++i) // Place upscale number of points for each key
		{
			if(is_diffusionSpace(grid.template get<PHI_SDF_ECS_FULL>(key), b_low, b_up))
			{
				sparse_grid.template insertFlush<PHI_SDF_ECS_SPARSE>(key_upscale) = grid.template get<PHI_SDF_ECS_FULL>(key);
			}
			++dom_upscale
		}


		++dom;
	}

}
#endif

template <
size_t PHI_SDF_ECS_FULL,
size_t PHI_SDF_SHELL_FULL,
size_t PHI_SDF_ECS_SPARSE,
size_t PHI_SDF_SHELL_SPARSE,
typename grid_type, 
typename sparse_in_type, 
typename T>
void get_diffusion_domain_sparse_grid_with_shell(grid_type & grid_ecs,
									  grid_type & grid_shell,
									  sparse_in_type & sparse_grid,
									  const T b_low,
									  const T b_up)
{
	auto dom = grid_ecs.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		if(is_diffusionSpace(grid_ecs.template get<PHI_SDF_ECS_FULL>(key), b_low, b_up))
		{
			sparse_grid.template insertFlush<PHI_SDF_ECS_SPARSE>(key)   = grid_ecs.template get<PHI_SDF_ECS_FULL>(key);
			sparse_grid.template insertFlush<PHI_SDF_SHELL_SPARSE>(key) = grid_shell.template get<PHI_SDF_SHELL_FULL>(key);
		}
		++dom;
	}
}


#if 0
// Init c0
	auto dom = g_sparse.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		if (g_sparse.template getProp<PHI_SDF>(key) < 4 * g_sparse.spacing(x))
		{
			g_sparse.template getProp<CONC_N>(key) =
			        (1 - 0.25 * g_sparse.template getProp<PHI_SDF>(key)) / g_sparse.spacing(x);
		}
		++dom;
	}
	if(v_cl.rank() == 0) std::cout << "Initialized grid with smoothed c0 at the boundary." << std::endl;
	g_sparse.write(path_output + "g_init", FORMAT_BINARY);
#endif




template <typename point_type, typename y_margin_type>
auto distance_from_margin(point_type & coord, y_margin_type y_margin)
{
	return(coord.get(1) - y_margin);
}

template <typename point_type, typename y_margin_type, typename width_type>
bool is_source(point_type & coord, y_margin_type y_margin, width_type width_source)
{
	return (coord.get(1) >= y_margin
			&& distance_from_margin(coord, y_margin) <= width_source);
}

template <typename T>
static bool is_inner_surface(const T & phi_sdf, const T & lower_bound)
{
	const T EPSILON = std::numeric_limits<T>::epsilon();
	const T _b_low = lower_bound + EPSILON;
	return (phi_sdf > _b_low);
}

template <size_t PHI_SDF, size_t K_SOURCE, size_t K_SINK, typename grid_type, typename y_margin_type, typename width_type, typename k_type>
void init_reactionTerms(grid_type & grid,
                        width_type width_source,
						y_margin_type y_margin,
                        k_type k_source,
                        k_type k_sink,
                        int no_membrane_points=4)
{
	/*
	 * Setting the reaction terms according to their distance from the margin of the blastoderm along the dorsal-ventral axis
	 *
	 * */
	
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		Point<grid_type::dims, y_margin_type> coord = grid.getPos(key);
		
		if(grid.template get<PHI_SDF>(key) < no_membrane_points * grid.spacing(0)) // If point lies close to the cell surface
		{
			if(is_source(coord, y_margin, width_source)) // If point lies within width_source distance from the margin, set k_source
			{
				grid.template insertFlush<K_SOURCE>(key) = k_source;
			}
			else
			{
				grid.template insertFlush<K_SOURCE>(key) = 0.0; // If membrane point not in source, emission is zero
			}
			grid.template insertFlush<K_SINK>(key) = k_sink; // Have sink at all membranes
			
		}
		
		else // For the nodes that are further away from the membrane, set both reaction terms to zero
		{
			grid.template insertFlush<K_SOURCE>(key) = 0.0;
			grid.template insertFlush<K_SINK>(key) = 0.0;
		}
		++dom;
	}
}

template <
size_t PHI_SDF_ECS_SPARSE, 
size_t PHI_SDF_SHELL_SPARSE, 
size_t K_SOURCE, 
size_t K_SINK, 
typename grid_type, 
typename y_margin_type, 
typename width_type, 
typename k_type,
typename boundary_type>
void init_reactionTerms_with_shell(grid_type & grid,
                        width_type width_source,
						y_margin_type y_margin,
                        k_type k_source,
                        k_type k_sink,
                        boundary_type b_low_shell,
                        int no_membrane_points=4)
{
	/*
	 * Setting the reaction terms according to their distance from the margin of the blastoderm along the dorsal-ventral axis
	 *
	 * */
	
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		Point<grid_type::dims, y_margin_type> coord = grid.getPos(key);
		
		if(
		grid.template get<PHI_SDF_ECS_SPARSE>(key) < no_membrane_points * grid.spacing(0) // If point is a surface point
		&& is_inner_surface(grid.template get<PHI_SDF_SHELL_SPARSE>(key), b_low_shell) // and this surface is not towards the yolk or EVL
		) 
		{
			if(is_source(coord, y_margin, width_source)) // If point lies within width_source distance from the margin, set k_source
			{
				grid.template insertFlush<K_SOURCE>(key) = k_source;
			}
			else
			{
				grid.template insertFlush<K_SOURCE>(key) = 0.0; // If membrane point not in source, emission is zero
			}
			grid.template insertFlush<K_SINK>(key) = k_sink; // Have sink at all membranes
			
		}
		
		else // For the nodes that are further away from the membrane or belong to the outer surface, set both reaction terms to zero
		{
			grid.template insertFlush<K_SOURCE>(key) = 0.0;
			grid.template insertFlush<K_SINK>(key) = 0.0;
		}
		++dom;
	}
}


template <size_t PHI_SDF, size_t K_SOURCE, size_t K_SINK, typename grid_type, typename coord_type, typename k_type>
void init_reactionTerms_smoothed(grid_type & grid,
                        coord_type y_start,
                        coord_type y_peak,
                        coord_type y_stop,
                        k_type k_source,
                        k_type k_sink,
						coord_type no_membrane_points=4,
						coord_type smoothing=0.25)
{
	/*
	 * Setting the reaction terms according to the position along the animal-vegetal axis = y-axis
	 *
	 * */
	const int x = 0;
	const int y = 1;
	const int z = 2;
	
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		Point<grid_type::dims, coord_type> coord = grid.getPos(key);
//		k_type sdf_based_smoothing = (1 - smoothing * grid.template getProp<PHI_SDF>(key)) / grid.spacing(x);
		k_type sdf_based_smoothing = 1.0;
		if(grid.template get<PHI_SDF>(key) < no_membrane_points * grid.spacing(0)) // If grid node lies close to the
			// membrane
		{
			// If membrane on animal side of the peak, emission strength increases with y
			if (coord[y] > y_start
			&& coord[y] < y_peak + std::numeric_limits<coord_type>::epsilon())
			{
				grid.template insertFlush<K_SOURCE>(key) = k_source * (coord[y] - y_start) * sdf_based_smoothing;
			}
			// Else if membrane on vegetal side of the peak, emission strength decreases with y as moving towards the
			// boundary to the yolk
			else if (coord[y] >= y_peak - std::numeric_limits<coord_type>::epsilon()
			&& coord[y] < y_stop - std::numeric_limits<coord_type>::epsilon())
			{
				grid.template insertFlush<K_SOURCE>(key) = k_source * (y_stop - coord[y]) * sdf_based_smoothing;
			}
			// Else source is zero
			else grid.template insertFlush<K_SOURCE>(key) = 0.0;
			grid.template insertFlush<K_SINK>(key) = k_sink; // Have sink at all membranes
		}
		
		else // For the nodes that are further away from the membrane, set the reaction terms to zero
		{
			grid.template insertFlush<K_SOURCE>(key) = 0.0;
			grid.template insertFlush<K_SINK>(key) = 0.0;
		}
		++dom;
	}
}



#endif //SPARSEGRID_DIFFUSIONSPACE_SPARSEGRID_HPP
