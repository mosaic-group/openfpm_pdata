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
#include "VTKWriter/VTKWriter.hpp"
#include "nn_processor.hpp"

/*! \brief structure that store and compute the internal and external local ghost box
 *
 * \tparam dim is the dimensionality of the physical domain we are going to decompose.
 * \tparam T type of the space we decompose, Real, Integer, Complex ...
 *
 * \see CartDecomposition
 *
 */
template<unsigned int dim, typename T>
class ie_loc_ghost
{
	//! It contain the calculated local ghost boxes
	openfpm::vector<lBox_dom<dim,T>> loc_ghost_box;

	//! temporal added sub-domains
	openfpm::vector<Box_loc_sub<dim,T>> sub_domains_tmp;

	/*! \brief Create the external local ghost boxes
	 *
	 * \param ghost part
	 * \param sub_domains list of local sub-domains
	 * \param sub_domains_prc list of sub-domains from the neighborhood processors
	 *
	 */
	void create_loc_ghost_ebox(Ghost<dim,T> & ghost,
			                   openfpm::vector<SpaceBox<dim,T>> & sub_domains,
							   openfpm::vector<Box_loc_sub<dim,T>> & sub_domains_prc)
	{
		comb<dim> zero;
		zero.zero();

		loc_ghost_box.resize(sub_domains.size());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			SpaceBox<dim,T> sub_with_ghost = sub_domains.get(i);

			// enlarge the sub-domain with the ghost
			sub_with_ghost.enlarge(ghost);

			// intersect with the other local sub-domains
			for (size_t j = 0 ; j < sub_domains_prc.size() ; j++)
			{
				size_t rj = sub_domains_prc.get(j).sub;

				if (rj == i && sub_domains_prc.get(j).cmb == zero)
					continue;

				::Box<dim,T> bi;

				bool intersect = sub_with_ghost.Intersect(sub_domains_prc.get(j).bx,bi);

				if (intersect == true)
				{
					Box_sub<dim,T> b;
					b.sub = rj;
					b.bx = bi;
					b.cmb = sub_domains_prc.get(j).cmb;

					// local external ghost box
					loc_ghost_box.get(i).ebx.add(b);

					// search this box in the internal box of the sub-domain j
					for (size_t k = 0; k < loc_ghost_box.get(rj).ibx.size() ; k++)
					{
						if (loc_ghost_box.get(rj).ibx.get(k).sub == i && loc_ghost_box.get(rj).ibx.get(k).cmb == sub_domains_prc.get(j).cmb.operator-())
						{
							loc_ghost_box.get(rj).ibx.get(k).k = loc_ghost_box.get(i).ebx.size()-1;
							break;
						}
					}
				}
			}
		}
	}

	/*! \brief Create the internal local ghost boxes
	 *
	 * \param ghost part
	 * \param sub_domains local sub-domains
	 * \param sub_domains_prc list of sub-domains from the neighborhood processors
	 *
	 */
	void create_loc_ghost_ibox(Ghost<dim,T> & ghost,
			                   openfpm::vector<SpaceBox<dim,T>> & sub_domains,
							   openfpm::vector<Box_loc_sub<dim,T>> & sub_domains_prc)
	{
		comb<dim> zero;
		zero.zero();

		loc_ghost_box.resize(sub_domains.size());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			// intersect with the others local sub-domains
			for (size_t j = 0 ; j < sub_domains_prc.size() ; j++)
			{
				SpaceBox<dim,T> sub_with_ghost = sub_domains_prc.get(j).bx;
				size_t rj = sub_domains_prc.get(j).sub;

				// Avoid to intersect the box with itself
				if (rj == i && sub_domains_prc.get(j).cmb == zero)
					continue;

				// enlarge the sub-domain with the ghost
				sub_with_ghost.enlarge(ghost);

				::Box<dim,T> bi;

				bool intersect = sub_with_ghost.Intersect(::SpaceBox<dim,T>(sub_domains.get(i)),bi);

				if (intersect == true)
				{
					Box_sub_k<dim,T> b;
					b.sub = rj;
					b.bx = bi;
					b.k = -1;
					b.cmb = sub_domains_prc.get(j).cmb;

					loc_ghost_box.get(i).ibx.add(b);
				}
			}
		}
	}

	/*! \brief In case of periodic boundary conditions we replicate the sub-domains at the border
	 *
	 * \param sub_domains list of sub-domains
	 * \param domain Domain box
	 * \param ghost part
	 * \param bc boundary conditions
	 *
	 */
	void applyBC(openfpm::vector<Box_loc_sub<dim,T>> & sub_domains, const Box<dim,T> & domain, const Ghost<dim,T> & ghost, const size_t (&bc)[dim])
	{
		HyperCube<dim> hyp;

		// first we create boxes at the border of the domain used to detect the sub-domain
		// that must be adjusted, each of this boxes define a shift in case of periodic boundary condition
		for (long int i = dim-1 ; i >= 0 ; i--)
		{
			std::vector<comb<dim>> cmbs = hyp.getCombinations_R(i);

			for (size_t j = 0 ; j < cmbs.size() ; j++)
			{
				if (nn_prcs<dim,T>::check_valid(cmbs[j],bc) == false)
					continue;

				Box<dim,T> bp;
				Point<dim,T> shift;

				for (size_t k = 0 ; k < dim ; k++)
				{
					switch (cmbs[j][k])
					{
					case 1:
						bp.setLow(k,domain.getHigh(k)+ghost.getLow(k));
						bp.setHigh(k,domain.getHigh(k));
						shift.get(k) = -domain.getHigh(k)+domain.getLow(k);
						break;
					case 0:
						bp.setLow(k,domain.getLow(k));
						bp.setHigh(k,domain.getHigh(k));
						shift.get(k) = 0;
						break;
					case -1:
						bp.setLow(k,domain.getLow(k));
						bp.setHigh(k,ghost.getHigh(k));
						shift.get(k) = domain.getHigh(k)-domain.getLow(k);
						break;
					}
				}

				// Detect all the sub-domain involved, shift them and add to the list
				// Detection is performed intersecting the sub-domains with the ghost
				// parts near the domain borders
				for (size_t k = 0 ; k < sub_domains.size() ; k++)
				{
					Box<dim,T> sub = sub_domains.get(k).bx;
					Box<dim,T> b_int;

					if (sub.Intersect(bp,b_int) == true)
					{
						sub += shift;
						add_subdomain(Box_loc_sub<dim,T>(sub,k,cmbs[j]));
					}
				}
			}
		}

		flush(sub_domains);
	}



	/*! \brief add sub-domains to a temporal list
	 *
	 * \param bx Box to add
	 *
	 */
	inline void add_subdomain(const Box_loc_sub<dim,T> & bx)
	{
		sub_domains_tmp.add(bx);
	}

	/*! \brief Flush the temporal added sub-domain to the sub-domain list
	 *
	 * \param sub_domains to add (In general they come from mirroring periodic
	 *        boundary conditions)
	 *
	 */
	void flush(openfpm::vector<Box_loc_sub<dim,T>> & sub_domains)
	{
		for (size_t i = 0 ; i < sub_domains_tmp.size() ; i++)
		{
			sub_domains.add(sub_domains_tmp.get(i));
		}

		sub_domains_tmp.clear();
	}

public:

	/*! \brief Create external and internal local ghosts
	 *
	 * \param sub_domains list of local sub-domains
	 * \param domain simulation domain
	 * \param ghost boundary
	 * \param bc Boundary conditions
	 *
	 */
	void create(openfpm::vector<SpaceBox<dim,T>> & sub_domains, Box<dim,T> & domain , Ghost<dim,T> & ghost , const size_t (&bc)[dim] )
	{
		// It will store local sub-domains + borders
		openfpm::vector<Box_loc_sub<dim,T>> sub_domains_prc;

		comb<dim> zero;
		zero.zero();

		// Copy sub_domains into sub_domains_prc
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			Box_loc_sub<dim,T> bls(SpaceBox<dim,T>(sub_domains.get(i)),i,zero);
			sub_domains_prc.add(bls);
			sub_domains_prc.last().sub = i;
		}

		applyBC(sub_domains_prc,domain,ghost,bc);

		create_loc_ghost_ibox(ghost,sub_domains,sub_domains_prc);
		create_loc_ghost_ebox(ghost,sub_domains,sub_domains_prc);
	}

	//! Default constructor
	ie_loc_ghost()	{};

	//! Constructor from another ie_loc_ghost
	ie_loc_ghost(const ie_loc_ghost<dim,T> & ilg)
	{
		this->operator=(ilg);
	};

	//! Constructor from temporal ie_loc_ghost
	ie_loc_ghost(ie_loc_ghost<dim,T> && ilg)
	{
		this->operator=(ilg);
	}

	/*! \brief copy the ie_loc_ghost
	 *
	 * \param ilg object to copy
	 *
	 * \return itself
	 *
	 */
	ie_loc_ghost<dim,T> & operator=(const ie_loc_ghost<dim,T> & ilg)
	{
		loc_ghost_box = ilg.loc_ghost_box;
		return *this;
	}

	/*! \brief copy the ie_loc_ghost
	 *
	 * \param ilg object to copy
	 *
	 * \return itself
	 *
	 */
	ie_loc_ghost<dim,T> & operator=(ie_loc_ghost<dim,T> && ilg)
	{
		loc_ghost_box.swap(ilg.loc_ghost_box);
		return *this;
	}

	/*! \brief Get the number of local sub-domains
	 *
	 * \return the number of local sub-domains
	 *
	 *
	 */
	inline size_t getNLocalSub()
	{
		return loc_ghost_box.size();
	}

	/*! \brief Get the number of external local ghost box for each sub-domain
	 *
	 * \param id sub-domain id
	 *
	 * \return the number of external ghost box
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
	 * \return the number of internal ghost box
	 *
	 */
	inline size_t getLocalNIGhost(size_t id)
	{
		return loc_ghost_box.get(id).ibx.size();
	}

	/*! \brief For the sub-domain i intersected with the sub-domain j enlarged, the associated
	 *       external ghost box is located in getLocalEGhostBox(j,k) with
	 *       getLocalIGhostSub(j,k) == i, this function return k
	 *
	 * \param i
	 * \param j
	 *
	 * \return k
	 *
	 */
	inline size_t getLocalIGhostE(size_t i, size_t j)
	{
		return loc_ghost_box.get(i).ibx.get(j).k;
	}

	/*! \brief Get the j internal local ghost box for the i sub-domain
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
		return loc_ghost_box.get(i).ibx.get(j).bx;
	}

	/*! \brief Get the j internal local ghost box boundary position for the i sub-domain of the local processor
	 *
	 * \note For the sub-domain i intersected with the sub-domain j enlarged, the associated
	 *       external ghost box is located in getLocalIGhostBox(j,k) with
	 *       getLocalIGhostSub(j,k) == i
	 *
	 * To get k use getLocalIGhostE
	 *
	 * \see getLocalIGhostE
	 *
	 * Some of the intersection boxes has special position, because they are at the boundary, this function
	 * return their position at the border
	 *
		\verbatim

															[1,1]
			+---------+------------------------+---------+
			| (1,-1)  |                        | (1,1)   |
			|   |     |    (1,0) --> 7         |   |     |
			|   v     |                        |   v     |
			|   6     |                        |   8     |
			+--------------------------------------------+
			|         |                        |         |
			|         |                        |         |
			|         |                        |         |
			| (-1,0)  |                        | (1,0)   |
			|    |    |                        |   |     |
			|    v    |      (0,0) --> 4       |   v     |
			|    3    |                        |   5     |
			|         |                        |         |
		 	|         |                        |         |
			|         |                        |         |
			|         |                        |         |
			|         |                        |         |
			|         |                        |         |
			+--------------------------------------------+
			| (-1,-1) |                        | (-1,1)  |
			|    |    |   (-1,0) --> 1         |    |    |
			|    v    |                        |    v    |
			|    0    |                        |    2    |
			+---------+------------------------+---------+


		\endverbatim
	 *
	 * \param i sub-domain
	 * \param j box
	 * \return the box
	 *
	 */
	inline const comb<dim> & getLocalIGhostPos(size_t i, size_t j) const
	{
		return loc_ghost_box.get(i).ibx.get(j).cmb;
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
		return loc_ghost_box.get(i).ebx.get(j).bx;
	}

	/*! \brief Get the j external local ghost box for the local processor
	 *
	 * \param i sub-domain
	 * \param j box
	 * \return the box
	 *
	 */
	inline const comb<dim> & getLocalEGhostPos(size_t i, size_t j) const
	{
		return loc_ghost_box.get(i).ebx.get(j).cmb;
	}

	/*! \brief Considering that sub-domain has N internal local ghost box identified
	 *         with the 0 <= k < N that come from the intersection of 2 sub-domains i and j
	 *         where j is enlarged, given the sub-domain i and the id k to identify the local internal ghost,
	 *          it return the id k of the other sub-domain that produced the intersection
	 *
	 * \param i sub-domain
	 * \param k id
	 * \return j
	 *
	 */
	inline size_t getLocalIGhostSub(size_t i, size_t k) const
	{
		return loc_ghost_box.get(i).ibx.get(k).sub;
	}

	/*! \brief Considering that sub-domain has N external local ghost box identified
	 *         with the 0 <= k < N that come from the intersection of 2 sub-domains i and j
	 *         where i is enlarged, given the sub-domain i and the id k of the external box,
	 *         it return the id of the other sub-domain that produced the intersection
	 *
	 * \param i sub-domain
	 * \param k id
	 * \return j
	 *
	 */
	inline size_t getLocalEGhostSub(size_t i, size_t k) const
	{
		return loc_ghost_box.get(i).ebx.get(k).sub;
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
	 * \return true if the file is written correctly
	 *
	 */
	bool write(std::string output, size_t p_id) const
	{
		// Copy the Box_sub_k into a vector of boxes
		openfpm::vector<openfpm::vector<Box<dim,T>>> vv5;

		for (size_t p = 0 ; p < loc_ghost_box.size() ; p++)
		{
			vv5.add();
			for (size_t i = 0 ; i < loc_ghost_box.get(p).ibx.size() ; i++)
				vv5.last().add(loc_ghost_box.get(p).ibx.get(i).bx);
		}

		//! local_internal_ghost_X.vtk internal local ghost boxes for the local processor (X)
		VTKWriter<openfpm::vector<Box<dim,T>>,VECTOR_BOX> vtk_box5;
		for (size_t p = 0 ; p < vv5.size() ; p++)
		{
			vtk_box5.add(vv5.get(p));
		}
		vtk_box5.write(output + std::string("local_internal_ghost_") + std::to_string(p_id) + std::string(".vtk"));

		// Copy the Box_sub_k into a vector of boxes
		openfpm::vector<openfpm::vector<Box<dim,T>>> vv6;

		for (size_t p = 0 ; p < loc_ghost_box.size() ; p++)
		{
			vv6.add();
			for (size_t i = 0 ; i < loc_ghost_box.get(p).ebx.size() ; i++)
				vv6.last().add(loc_ghost_box.get(p).ebx.get(i).bx);
		}

		//! local_external_ghost_X.vtk external local ghost boxes for the local processor (X)
		VTKWriter<openfpm::vector<Box<dim,T>>,VECTOR_BOX> vtk_box6;
		for (size_t p = 0 ; p < vv6.size() ; p++)
		{
			vtk_box6.add(vv6.get(p));
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
				{
					std::cout << "No ibx link" << "\n";
					return false;
				}
			}
		}

		return true;
	}

	/*! \brief Check if the ie_loc_ghosts contain the same information
	 *
	 * \param ilg Element to check
	 *
	 * \return true if they match
	 *
	 */
	bool is_equal(ie_loc_ghost<dim,T> & ilg)
	{
		if (ilg.loc_ghost_box.size() != loc_ghost_box.size())
			return false;

		// Explore all the subdomains
		for (size_t i = 0 ; i < loc_ghost_box.size() ; i++)
		{
			if (getLocalNIGhost(i) != ilg.getLocalNIGhost(i))
				return false;

			if (getLocalNEGhost(i) != ilg.getLocalNEGhost(i))
				return false;

			for (size_t j = 0 ; j < getLocalNIGhost(i) ; j++)
			{
				if (getLocalIGhostE(i,j) != ilg.getLocalIGhostE(i,j))
					return false;
				if (getLocalIGhostBox(i,j) != ilg.getLocalIGhostBox(i,j))
					return false;
				if (getLocalIGhostSub(i,j) != ilg.getLocalIGhostSub(i,j))
					return false;
			}
			for (size_t j = 0 ; j < getLocalNEGhost(i) ; j++)
			{
				if (getLocalEGhostBox(i,j) != ilg.getLocalEGhostBox(i,j))
					return false;
				if (getLocalEGhostSub(i,j) != ilg.getLocalEGhostSub(i,j))
					return false;
			}

		}

		return true;
	}



	/*! \brief Check if the ie_loc_ghosts contain the same information
	 * with the exception of the ghost part
	 *
	 * \param ilg Element to check
	 *
	 * \return true if the two objects are equal with the exception of the
	 *         ghost part
	 *
	 */
	bool is_equal_ng(ie_loc_ghost<dim,T> & ilg)
	{
		Box<dim,T> bt;

		if (ilg.loc_ghost_box.size() != loc_ghost_box.size())
			return false;

		// Explore all the subdomains
		for (size_t i = 0 ; i < loc_ghost_box.size() ; i++)
		{
			if (getLocalNIGhost(i) != ilg.getLocalNIGhost(i))
				return false;

			if (getLocalNEGhost(i) != ilg.getLocalNEGhost(i))
				return false;

			for (size_t j = 0 ; j < getLocalNIGhost(i) ; j++)
			{
				if (getLocalIGhostE(i,j) != ilg.getLocalIGhostE(i,j))
					return false;
				if (getLocalIGhostBox(i,j).Intersect(ilg.getLocalIGhostBox(i,j),bt) == false)
					return false;
				if (getLocalIGhostSub(i,j) != ilg.getLocalIGhostSub(i,j))
					return false;
			}
			for (size_t j = 0 ; j < getLocalNEGhost(i) ; j++)
			{
				if (getLocalEGhostBox(i,j).Intersect(ilg.getLocalEGhostBox(i,j),bt) == false)
					return false;
				if (getLocalEGhostSub(i,j) != ilg.getLocalEGhostSub(i,j))
					return false;
			}

		}

		return true;
	}

	/*! \brief Reset the ie_loc_ghost
	 *
	 */
	void reset()
	{
		loc_ghost_box.clear();
		sub_domains_tmp.clear();
	}
};


#endif /* SRC_DECOMPOSITION_GHOST_DEC_IE_GHOST_HPP_ */
