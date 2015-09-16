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
	inline size_t getSub() const
	{
		return g_c;
	}

	/*! \brief Get the key
	 *
	 * \return the local key
	 *
	 */
	inline grid_key_dx<dim> getKey() const
	{
		return key;
	}

	/*! \brief Create a new key moving the old one
	 *
	 * \param i dimension id
	 * \param s number of steps
	 *
	 * \return new key
	 *
	 */
	inline grid_dist_key_dx<dim> move(size_t i,size_t s)
	{
		grid_key_dx<dim> key = getKey();
		key.set_d(i,key.get(i) + s);
		return grid_dist_key_dx<dim>(getSub(),key);
	}

	inline grid_dist_key_dx(int g_c, grid_key_dx<dim> key)
	:g_c(g_c),key(key)
	{
	}
};

#endif
