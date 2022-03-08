#include "Space/SpaceBox.hpp"

#ifndef DECOMPOSITION_HPP_
#define DECOMPOSITION_HPP_

/**
 *
 * \brief class that store Internal part external and border part of a dataset
 *
 */

template <typename T> class dataDiv
{
	//! Border part of the data
	boost::shared_ptr<T> bord;
	//! internal part of your data
	boost::shared_ptr<T> inte;
	//! external part of your data
	boost::shared_ptr<T> ext;
};

/*! \brief This class define the domain decomposition interface
 *
 * This class define the domain decomposition interface, its main functionality
 * is to divide a domain in several subspaces
 *
 * \tparam T structure that store the dataset
 * \tparam S type of space is decomposing Real integer complex ...
 *
 */

template<typename T, typename S>
class Decomposition
{
	/*! \brief The the internal part of the data set, or the data that
	 * does not depend from the ghosts layers
	 *
	 * \return The internal part of the dataset
	 *
	 */
	virtual T getInternal();


	/*! Get the ghost part of the dataset
	 *
	 * \return The internal part of the dataset
	 *
	 */
	virtual T getBorder();

	/*! Get the external part of the dataset (outside the ghost)
	 *
	 * \return The external part of the dataset
	 *
	 */
	virtual T getExternal();

	//! divide the dataset from internal part and border
	virtual dataDiv<T> divide();

	//! Get the number of hyper-cube the space id is divided into
	virtual size_t getNHyperCube(size_t id);

	//! Get the hyper-cube margins
	virtual std::vector<T> & getHyperCube(size_t id, size_t id_c);

	//! destructor
	virtual ~Decomposition(){}
};

#endif
