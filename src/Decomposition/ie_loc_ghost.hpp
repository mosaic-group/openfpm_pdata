/*
 * ie_ghost.hpp
 *
 *  Created on: Aug 8, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_GHOST_DEC_IE_GHOST_HPP_
#define SRC_DECOMPOSITION_GHOST_DEC_IE_GHOST_HPP_

#include "Space/Shape/Box.hpp"
#include "Space/Ghost.hpp"
#include "Space/SpaceBox.hpp"
#include "common.hpp"
#include "VTKWriter.hpp"

/*! \brief structure that store and compute the internal and external local ghost box
 *
 * \see CartDecomposition
 *
 */

template<unsigned int dim, typename T>
class ie_loc_ghost
{
	openfpm::vector<lBox_dom<dim,T>> loc_ghost_box;

	// Save the ghost boundaries
//	Ghost<dim,T> ghost;

protected:

	/*! \brief Create the external local ghost boxes
	 *
	 * \param ghost margin to enlarge
	 * \param local sub-domain
	 *
	 */
	void create_loc_ghost_ebox(Ghost<dim,T> & ghost, openfpm::vector<SpaceBox<dim,T>> & sub_domains)
	{
		// Save the ghost
//		this->ghost = ghost;

		loc_ghost_box.resize(sub_domains.size());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			SpaceBox<dim,T> sub_with_ghost = sub_domains.get(i);

			// enlarge the sub-domain with the ghost
			sub_with_ghost.enlarge(ghost);

			// intersect with the other local sub-domains
			for (size_t j = 0 ; j < sub_domains.size() ; j++)
			{
				if (i == j)
					continue;

				::Box<dim,T> bi;

				bool intersect = sub_with_ghost.Intersect(::SpaceBox<dim,T>(sub_domains.get(j)),bi);

				if (intersect == true)
				{
					Box_sub<dim,T> b;
					b.sub = j;
					b = bi;

					// local external ghost box
					loc_ghost_box.get(i).ebx.add(b);

					// search this box in the internal box of the sub-domain j
					for (size_t k = 0; k < loc_ghost_box.get(j).ibx.size() ; k++)
					{
						if (loc_ghost_box.get(j).ibx.get(k).sub == i)
						{
							loc_ghost_box.get(j).ibx.get(k).k = loc_ghost_box.get(i).ebx.size()-1;
							break;
						}
					}
				}
			}
		}
	}

	/*! \brief Create the internal local ghost boxes
	 *
	 * \param ghost margin to enlarge
	 *
	 */
	void create_loc_ghost_ibox(Ghost<dim,T> & ghost, openfpm::vector<SpaceBox<dim,T>> & sub_domains)
	{
		loc_ghost_box.resize(sub_domains.size());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			// intersect with the others local sub-domains
			for (size_t j = 0 ; j < sub_domains.size() ; j++)
			{
				if (i == j)
					continue;

				SpaceBox<dim,T> sub_with_ghost = sub_domains.get(j);
				// enlarge the sub-domain with the ghost
				sub_with_ghost.enlarge(ghost);

				::Box<dim,T> bi;

				bool intersect = sub_with_ghost.Intersect(::SpaceBox<dim,T>(sub_domains.get(i)),bi);

				if (intersect == true)
				{
					Box_sub_k<dim,T> b;
					b.sub = j;
					b = bi;
					b.k = -1;

					loc_ghost_box.get(i).ibx.add(b);
				}
			}
		}
	}

public:

	/*! \brief Get the number of external local ghost box for each sub-domain
	 *
	 * \param id sub-domain id
	 *
	 * \return the number of internal ghost box
	 *
	 */
	inline size_t getLocalNEGhost(size_t id)
	{
		return loc_ghost_box.get(id).ebx.size();
	}

	/*! \brief Get the number of internal local ghost box for each sub-domain
	 *
	 * \param id sub-domain id
	 *
	 * \return the number of external ghost box
	 *
	 */
	inline size_t getLocalNIGhost(size_t id)
	{
		return loc_ghost_box.get(id).ibx.size();
	}

	/*! \brief For the sub-domain i intersected with the sub-domain j enlarged, the associated
	 *       external ghost box is located in getLocalIGhostBox(j,k) with
	 *       getLocalIGhostSub(j,k) == i, this function return k
	 *
	 *
	 */
	inline size_t getLocalIGhostE(size_t i, size_t j)
	{
		return loc_ghost_box.get(i).ibx.get(j).k;
	}

	/*! \brief Get the j internal local ghost box for the i sub-domain of the local processor
	 *
	 * \note For the sub-domain i intersected with the sub-domain j enlarged, the associated
	 *       external ghost box is located in getLocalIGhostBox(j,k) with
	 *       getLocalIGhostSub(j,k) == i
	 *
	 * To get k use getLocalIGhostE
	 *
	 * \see getLocalIGhostE
	 *
	 * \param i sub-domain
	 * \param j box
	 * \return the box
	 *
	 */
	inline const ::Box<dim,T> & getLocalIGhostBox(size_t i, size_t j) const
	{
		return loc_ghost_box.get(i).ibx.get(j);
	}

	/*! \brief Get the j external local ghost box for the local processor
	 *
	 * \param i sub-domain
	 * \param j box
	 * \return the box
	 *
	 */
	inline const ::Box<dim,T> & getLocalEGhostBox(size_t i, size_t j) const
	{
		return loc_ghost_box.get(i).ebx.get(j);
	}

	/*! \brief Considering that sub-domain has N internal local ghost box identified
	 *         with the 0 <= k < N that come from the intersection of 2 sub-domains i and j
	 *         where j is enlarged, given the sub-domain i and the id k, it return the id of
	 *         the other sub-domain that produced the intersection
	 *
	 * \param i sub-domain
	 * \param k id
	 * \return the box
	 *
	 */
	inline size_t getLocalIGhostSub(size_t i, size_t j) const
	{
		return loc_ghost_box.get(i).ibx.get(j).sub;
	}

	/*! \brief Considering that sub-domain has N external local ghost box identified
	 *         with the 0 <= k < N that come from the intersection of 2 sub-domains i and j
	 *         where j is enlarged, given the sub-domain i and the id k, it return the id of
	 *         the other sub-domain that produced the intersection
	 *
	 * \param i sub-domain
	 * \param k id
	 * \return the box
	 *
	 */
	inline size_t getLocalEGhostSub(size_t i, size_t j) const
	{
		return loc_ghost_box.get(i).ebx.get(j).sub;
	}

	/*! \brief Write the decomposition as VTK file
	 *
	 * The function generate several files
	 *
	 * 5) local_internal_ghost_X.vtk internal local ghost boxes for the local processor (X)
	 * 6) local_external_ghost_X.vtk external local ghost boxes for the local processor (X)
	 *
	 * where X is the local processor rank
	 *
	 * \param output directory where to write the files
	 * \param p_id id of the local processor
	 *
	 */
	bool write(std::string output, size_t p_id) const
	{
		//! local_internal_ghost_X.vtk internal local ghost boxes for the local processor (X)
		VTKWriter<openfpm::vector_std<Box_sub_k<dim,T>>,VECTOR_BOX> vtk_box5;
		for (size_t p = 0 ; p < loc_ghost_box.size() ; p++)
		{
			vtk_box5.add(loc_ghost_box.get(p).ibx);
		}
		vtk_box5.write(output + std::string("local_internal_ghost_") + std::to_string(p_id) + std::string(".vtk"));

		//! local_external_ghost_X.vtk external local ghost boxes for the local processor (X)
		VTKWriter<openfpm::vector_std<Box_sub<dim,T>>,VECTOR_BOX> vtk_box6;
		for (size_t p = 0 ; p < loc_ghost_box.size() ; p++)
		{
			vtk_box6.add(loc_ghost_box.get(p).ebx);
		}
		vtk_box6.write(output + std::string("local_external_ghost_") + std::to_string(p_id) + std::string(".vtk"));

		return true;
	}

	/*! \brief function to check the consistency of the information of the decomposition
	 *
	 * \param n_sub Number of sub_domain
	 *
	 * \return false if is inconsistent
	 *
	 */
	bool check_consistency(size_t n_sub)
	{
		//! for each sub-domain
		for (size_t i = 0 ; i < loc_ghost_box.size() ; i++)
		{
			for (size_t j = 0 ; j < loc_ghost_box.get(i).ibx.size() ; j++)
			{
				if (loc_ghost_box.get(i).ibx.get(j).k == -1)
					return false;
			}
		}

		if (n_sub > 1)
		{
			for (size_t i = 0 ; i < loc_ghost_box.size() ; i++)
			{
				if (loc_ghost_box.get(i).ibx.size() == 0)
					return false;
				if (loc_ghost_box.get(i).ebx.size() == 0)
					return false;
			}
		}

		return true;
	}
};


#endif /* SRC_DECOMPOSITION_GHOST_DEC_IE_GHOST_HPP_ */
