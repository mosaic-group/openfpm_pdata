/*
 * basic.hpp
 *
 *  Created on: Aug 8, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_COMMON_HPP_
#define SRC_DECOMPOSITION_COMMON_HPP_

#define UNIQUE 1
#define MULTIPLE 2

#include "Vector/map_vector.hpp"

/*! \brief for each sub-domain box sub contain the real the sub-domain id
 *
 * When we apply boundary conditions real sub-domains are copied along the border
 * sub, contain the id of the real sub_domain
 *
 * \tparam dim Dimensionality of the box
 * \tparam T in witch space this box live
 *
 */
template<unsigned int dim, typename T>
struct Box_loc_sub
{
	//! Box defining the sub-domain (copied)
	Box<dim,T> bx;

	//! The id of the real domain
	size_t sub;

	//! in witch sector this sub-domain live
	comb<dim> cmb;

	//! Constructor
	Box_loc_sub()
	:sub(0)
	{
		cmb.zero();
	};

	//! Constructor from box, domain id and sector where it live
	Box_loc_sub(const Box<dim,T> & bx, size_t sub, const comb<dim> & cmb)
	:bx(bx),sub(sub),cmb(cmb)
	{};

	//! Set the sub-domain box coordinates
	Box_loc_sub operator=(const Box<dim,T> & box)
	{
		::Box<dim,T>::operator=(box);

		return *this;
	}


};

/*! It contain a box definition and from witch sub-domain it come from (in the local processor list)
 * and an unique if across adjacent processors (for communication)
 *
 * If the box come from the intersection of an expanded sub-domain and a sub-domain
 *
 * Assuming we are considering the near processors i (0 to getNNProcessors())
 *
 *
 */
template<unsigned int dim, typename T>
struct Box_sub
{
	//! Internal ghost box definition
	Box<dim,T> bx;

	//! Domain id
	size_t sub;

	//! see ebx_ibx_form in ie_ghost for the meaning
	size_t id;

	//! see getNearSubdomainsRealId in nn_prcs
	size_t r_sub;

	//! see ie_ghost follow sector explanation
	comb<dim> cmb;

	//! Constructor reset cmb
	Box_sub()
	{
		cmb.zero();
	}
};

//! Particular case for local internal ghost boxes
template<unsigned int dim, typename T>
struct Box_sub_k
{
	//! extension of this local internal ghost box
	Box<dim,T> bx;

	//! Domain id
	size_t sub;

	//! Where this sub_domain live
	comb<dim> cmb;

	//! k \see getLocalGhostIBoxE
	long int k;

	Box_sub_k()
	:sub(0),k(-1)
	{
		cmb.zero();
	}
};

template<unsigned int dim,typename T>
struct Box_map
{
	Box<dim,T> box;

	long int prc;

	static bool noPointers()
	{
		return true;
	}
};

//! Case for local ghost box
template<unsigned int dim, typename T>
struct lBox_dom
{
	//! Intersection between the local sub-domain enlarged by the ghost and the contiguous processor
	//! sub-domains (External ghost)
	openfpm::vector_std< Box_sub_k<dim,T> > ebx;

	//! Intersection between the contiguous processor sub-domain enlarged by the ghost with the
	//! local sub-domain (Internal ghost)
	openfpm::vector_std< Box_sub_k<dim,T>> ibx;
};

//! Case for local external ghost box
template<unsigned int dim, typename T>
struct Box_proc
{
	//! Intersection between the local sub-domain enlarged by the ghost and the contiguous processor
	//! sub-domains (External ghost)
	openfpm::vector<::Box<dim,T>> bx;

	//! Intersection between the contiguous processor sub-domain enlarged by the ghost with the
	//! local sub-domain (Internal ghost)
	openfpm::vector<::Box<dim,T>> nbx;


	//! processor
	size_t proc;
};

//! Case for external ghost box
template<unsigned int dim, typename T>
struct Box_dom
{
	//! Intersection between the local sub-domain enlarged by the ghost and the contiguous processor
	//! sub-domains (External ghost)
	openfpm::vector_std< Box_sub<dim,T> > ebx;

	//! Intersection between the contiguous processor sub-domain enlarged by the ghost with the
	//! local sub-domain (Internal ghost)
	openfpm::vector_std< Box_sub<dim,T> > ibx;
};

// It store the sub-domain sent by the near processors
template<unsigned int dim, typename T>
struct N_box
{
	//! id of the processor in the nn_processor list (local processor id)
	size_t id;

	//! near processor sub-domains
	typename openfpm::vector<::Box<dim,T>> bx;

	//! near processor sector position (or where they live outside the domain)
	openfpm::vector<comb<dim>> pos;

	//! Number of real sub-domains or sub-domain in the central sector
	size_t n_real_sub;

	//! When a sub-domain is not in the central sector, it mean that has been created
	//! because of periodicity in a non central sector. Any sub-domain not in the central
	//! sector is linked to one sub-domain in the central sector
	openfpm::vector<size_t> r_sub;

	//! Default constructor
	N_box()
	:id((size_t)-1),n_real_sub(0)
	{};

	//! Copy constructor
	N_box(const N_box<dim,T> & b)
	:id((size_t)-1),n_real_sub(0)
	{
		this->operator=(b);
	}

	//! Copy constructor
	N_box(N_box<dim,T> && b)
	:id((size_t)-1),n_real_sub(0)
	{
		this->operator=(b);
	}

	/*! \brief Copy the element
	 *
	 * \param ele element to copy
	 *
	 * \return itself
	 *
	 */
	N_box<dim,T> & operator=(const N_box<dim,T> & ele)
	{
		id = ele.id;
		bx = ele.bx;
		pos = ele.pos;
		n_real_sub = ele.n_real_sub;
		r_sub = ele.r_sub;

		return * this;
	}

	/*! \brief Copy the element
	 *
	 * \param ele element to copy
	 *
	 * \return itself
	 *
	 */
	N_box<dim,T> & operator=(N_box<dim,T> && ele)
	{
		id = ele.id;
		bx.swap(ele.bx);
		pos.swap(ele.pos);
		n_real_sub = ele.n_real_sub;
		r_sub.swap(ele.r_sub);

		return * this;
	}

	/*! \brief Compare two N_box object
	 *
	 * \param ele element to compare with
	 *
	 * \return true if they match
	 *
	 */
	bool operator==(const N_box<dim,T> & ele) const
	{
		if (id != ele.id)
			return false;

		if (pos != ele.pos)
			return false;

		if (r_sub != ele.r_sub)
			return false;

		if (n_real_sub != ele.n_real_sub)
			return false;

		return bx == ele.bx;
	}

	/*! \brief Compare two N_box object
	 *
	 * \param ele element to compare with
	 *
	 * \return true if they match
	 *
	 */
	bool operator!=(const N_box<dim,T> & ele) const
	{
		return ! this->operator==(ele);
	}
};

//! It store all the boxes of the near processors in a linear array
template<unsigned int dim, typename T>
struct p_box
{
	//! Box that identify the intersection of the ghost of the near processor with the
	//! processor sub-domain
	::Box<dim,T> box;
	//! local processor id
	size_t lc_proc;
	//! processor rank
	size_t proc;

	//! shift vector id
	size_t shift_id;

	/*! \brief Check if two p_box are the same
	 *
	 * \param pb box to check
	 *
	 * \return true if they match
	 *
	 */
	bool operator==(const p_box & pb)
	{
		return pb.lc_proc == lc_proc;
	}
};

#endif /* SRC_DECOMPOSITION_COMMON_HPP_ */
