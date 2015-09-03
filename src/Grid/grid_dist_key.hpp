#ifndef GRID_DIST_KEY_DX_HPP
#define GRID_DIST_KEY_DX_HPP

/*! \brief Grid key for a distributed grid
 *
 * It contain from which local sub-domain grid come from, and the local grid_key_dx
 *
 */

template<unsigned int dim>
class grid_dist_key_dx
{
	//! grid list counter

	size_t g_c;

	//! Local grid iterator

	grid_key_dx<dim> key;

public:

	/*! \brief Get the local grid
	 *
	 * \return the id of the local grid
	 *
	 */
	size_t getSub() const
	{
		return g_c;
	}

	/*! \brief Get the key
	 *
	 * \return the local key
	 *
	 */
	grid_key_dx<dim> getKey() const
	{
		return key;
	}

	grid_dist_key_dx(int g_c, grid_key_dx<dim> key)
	:g_c(g_c),key(key)
	{
	}

	std::string to_string()
	{
		std::stringstream str;

		str << "sub_domain=" << g_c << " ";

		for (size_t i = 0 ; i < dim ; i++)
			str << "x[" << i << "]=" << key.get(i) << " ";

		str << "\n";

		return str.str();
	}
};

#endif
