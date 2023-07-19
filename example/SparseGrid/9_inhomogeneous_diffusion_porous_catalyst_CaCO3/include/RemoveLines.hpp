//
// Created by jstark on 07.10.21.
//

#ifndef SUSSMAN_REDISTANCING_REMOVELINES_HPP
#define SUSSMAN_REDISTANCING_REMOVELINES_HPP

/**@brief Checks for thin lines and removes them by setting the value of the indicator function accordingly.
 *
 * @tparam grid_type Template type of OpenFPM grid.
 * @tparam Phi_0_grid Property that contains the indicator or level-set function.
 * @tparam sign_inside Sign that indicates that a node lies inside the object. Set to +1 or -1 depending on the
 * convention you use.
 * @param grid Grid that stores the indicator or level-set function.
 */
template <size_t Phi_0_grid, int sign_inside, typename grid_type>
void removeLines(grid_type & grid)
{
	if(sign_inside == 0)
	{
		std::cout << "sign_inside was set to 0. However, it must be set to +1 or -1 depending on the level-set "
					 "convention used. Aborting..." << std::endl;
		abort();
	}
	for(int repeat = 0; repeat < 2; ++repeat)
	{
		grid.template ghost_get<Phi_0_grid>(KEEP_PROPERTIES);
		if (sign_inside > 0)
		{
			auto dom = grid.getDomainIterator();
			while(dom.isNext())
			{
				auto key = dom.get();
				if (grid.template getProp<Phi_0_grid>(key) > 0)
				{
					for (int d = 0; d < grid_type::dims; ++d)
					{
						if (grid.template getProp<Phi_0_grid>(key.move(d, -1)) < 0
								&& grid.template getProp<Phi_0_grid>(key.move(d, +1)) < 0)
						{
							grid.template getProp<Phi_0_grid>(key) = grid.template getProp<Phi_0_grid>(key.move(d, +1));
							break;
						}
					}
				}
				++dom;
			}
		}
		if (sign_inside < 0)
		{
			auto dom = grid.getDomainIterator();
			while(dom.isNext())
			{
				auto key = dom.get();
				if (grid.template getProp<Phi_0_grid>(key) < 0)
				{
					for (int d = 0; d < grid_type::dims; ++d)
					{
						if (grid.template getProp<Phi_0_grid>(key.move(d, -1)) > 0
								&& grid.template getProp<Phi_0_grid>(key.move(d, +1)) > 0)
						{
							grid.template getProp<Phi_0_grid>(key) = grid.template getProp<Phi_0_grid>(key.move(d, +1));
							break;
						}
					}
				}
				++dom;
			}
		}
		
	}
}

template <typename T>
bool is_inside(T phi)
{
	return phi >= 0 - std::numeric_limits<T>::epsilon();
}

template <typename T>
bool is_outside(T phi)
{
	return phi < 0 + std::numeric_limits<T>::epsilon();
}

/**@brief Checks for thin lines and thin spaces in between surfaces of 1dx thickness and removes them by setting the
 * value of the indicator function equal accordingly.
 *
 * @tparam grid_type Template type of OpenFPM grid.
 * @tparam Phi_0_grid Property that contains the indicator or level-set function.
 * @param grid Grid that stores the indicator or level-set function.
 */
template <size_t Phi_0_grid, typename grid_type>
void removeLinesAndThinSpaces(grid_type & grid, int repeat=3)
{
	for(int k = 0; k < repeat; ++k)
	{
		grid.template ghost_get<Phi_0_grid>(KEEP_PROPERTIES);
		auto dom = grid.getDomainIterator();
		while(dom.isNext())
		{
			auto key = dom.get();
			for (int d = 0; d < grid_type::dims; ++d)
			{
				if (is_inside(grid.template getProp<Phi_0_grid>(key))) // If inside but neighbors are outside
				{
					if (is_outside(grid.template getProp<Phi_0_grid>(key.move(d, -1)))
						&& is_outside(grid.template getProp<Phi_0_grid>(key.move(d, +1))))
					{
						grid.template getProp<Phi_0_grid>(key) = grid.template getProp<Phi_0_grid>(key.move(d, +1));
						break;
					}
				}
				else if (is_outside(grid.template getProp<Phi_0_grid>(key))) // If outside but neighbors are inside
				{
					if (is_inside(grid.template getProp<Phi_0_grid>(key.move(d, -1)))
						&& is_inside(grid.template getProp<Phi_0_grid>(key.move(d, +1))))
					{
						grid.template getProp<Phi_0_grid>(key) = grid.template getProp<Phi_0_grid>(key.move(d, +1));
						break;
					}
				}
			}
			++dom;
		}
	}
}
#endif //SUSSMAN_REDISTANCING_REMOVELINES_HPP
