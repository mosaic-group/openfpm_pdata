/*
 * CartDecomposition.cpp
 *
 *  Created on: Aug 15, 2014
 *      Author: Pietro Incardona
 */

#include "CartDecomposition.hpp"



/*! \brief The the bulk part of the data set, or the data that does not depend
 *  from the ghosts layers
 *
 * The the bulk part of the data set, or the data that does not depend from the
 *  ghosts layers
 *
 */

/*template<typename T> T CartDecomposition<T>::getBulk(T data)
{
	// for each element in data

	for (size_t i = 0; i < data.size() ; i++)
	{
		if (localSpace.isInside())
	}

}

template<typename T> T CartDecomposition<T>::getInternal()
{

}*/

/*! \brief Check if is border or bulk
 *
 * \param neighboorhood define the neighboorhood of all the points
 * \return true if border, false if bulk
 *
 */

bool borderOrBulk(neighborhood & nb)
{
	device::grid<1,size_t> nbr = nb.next();

	// check the neighborhood

	// get neighborhood iterator

	grid_key_dx_iterator<dim> iterator_nbr = nbr.getIterator();

	while (iterator_nbr.hasNext())
	{
		grid_key_dx key_nbr = iterator_nbr.next();

		// check if the neighboorhood is internal

		if(subspace.isBound(data.template get<Point::x>(key_nbr)) == false)
		{
			// it is border

			return true;

			ret.bord.push_back(key);
			break;
		}
	}

	return false;
}

/*! \brief This function divide the data set into bulk, border, external and internal part
 *
 * \tparam dim dimensionality of the structure storing your data
 *         (example if they are in 3D grid, has to be 3)
 * \tparam T type of object we are dividing
 * \tparam device type of layout selected
 * \param data 1-dimensional grid of point
 * \param nb define the neighborhood of all the points
 * \return a structure with the set of objects divided
 *
 */

template<unsigned int dim, typename T, template<typename> class layout, typename Memory, template<unsigned int, typename> class Domain, template<typename, typename, typename> class data_s>
dataDiv<T> CartDecomposition<dim,T,layout>::divide(device::grid<1,Point<dim,T>> & data, neighborhood & nb)
{
	//! allocate the 3 subset

	dataDiv<T> ret;

	ret.bord = new boost::shared_ptr<T>(new T());
	ret.inte = new boost::shared_ptr<T>(new T());
	ret.ext = new boost::shared_ptr<T>(new T());

	//! get grid iterator

	grid_key_dx_iterator<dim> iterator = data.getIterator();

	//! we iterate trough all the set of objects

	while (iterator.hasNext())
	{
		grid_key_dx<dim> key = iterator.next();

		//! Check if the object is inside the subspace

		if (subspace.isBound(data.template get<Point<3,T>::x>(key)))
		{
			//! Check if the neighborhood is inside the subspace

			if (borderOrBulk(nb) == true)
			{
				// It is border

				ret.bord.push_back(key);
			}
			else
			{
				// It is bulk

				ret.bulk.push_back(key);
			}
		}
		else
		{
			//! it is external

			ret.ext.push_back(key);
		}
	}
}
