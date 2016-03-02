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
	Box<dim,T> bx;

	// Domain id
	size_t sub;

	// in witch sector this sub-domain live, when
	comb<dim> cmb;

	Box_loc_sub()
	{
		cmb.zero();
	};

	Box_loc_sub(const Box<dim,T> & bx, size_t sub, const comb<dim> & cmb)
	:bx(bx),sub(sub),cmb(cmb)
	{};

	template <typename Memory> Box_loc_sub(const Box_loc_sub<dim,T> & bls)
	{
		bx = bls.bx;
		this->sub = bls.sub;
	};

	Box_loc_sub operator=(const Box<dim,T> & box)
	{
		::Box<dim,T>::operator=(box);

		return *this;
	}


};

/*! It contain a box definition and from witch sub-domain it come from (in the local processor)
 * and an unique across adjacent processors (for communication)
 *
 * If the box come from the intersection of an expanded sub-domain and a sub-domain
 *
 * Assuming we are considering the near processors i (0 to getNNProcessors())
 *
 * ### external ghost box
 *
 * id = id_exp * N_non_exp + id_non_exp
 *
 * id_exp = the id in the vector proc_adj_box.get(i) of the expanded sub-domain (sent local sub-domains)
 *
 * id_non_exp = the id in the vector nn_processor_subdomains[i] of the sub-domain (received sub-domains from near processors)
 *
 * ### internal ghost box
 *
 * id = id_exp * N_non_exp + id_non_exp
 *
 * id_exp = the id in the vector nn_processor_subdomains[i] of the expanded sub-domain
 *
 * id_non_exp = the id in the vector proc_adj_box.get(i) of the sub-domain
 *
 */
template<unsigned int dim, typename T>
struct Box_sub
{
	Box<dim,T> bx;

	// Domain id
	size_t sub;

	// Id
	size_t id;

	Box_sub operator=(const Box<dim,T> & box)
	{
		bx = box;

		return *this;
	}


};

//! Particular case for local internal ghost boxes
template<unsigned int dim, typename T>
struct Box_sub_k
{
	Box<dim,T> bx;

	// Domain id
	size_t sub;

	// Where this sub_domain live
	comb<dim> cmb;

	//! k \see getLocalGhostIBoxE
	long int k;

	Box_sub_k()
	:k(-1)
	{
		cmb.zero();
	}

	Box_sub_k operator=(const Box<dim,T> & box)
	{
		bx = box;

		return *this;
	}

	// encap interface to make compatible with OpenFPM_IO
	template <int i> auto get() -> decltype( std::declval<Box<dim,T> *>()->template get<i>())
	{
		return ::Box<dim,T>::template get<i>();
	}
};

//! Case for local ghost box
template<unsigned int dim, typename T>
struct lBox_dom
{
	// Intersection between the local sub-domain enlarged by the ghost and the contiguous processor
	// sub-domains (External ghost)
	openfpm::vector_std< Box_sub<dim,T> > ebx;

	// Intersection between the contiguous processor sub-domain enlarged by the ghost with the
	// local sub-domain (Internal ghost)
	openfpm::vector_std< Box_sub_k<dim,T>> ibx;
};

template<unsigned int dim, typename T>
struct Box_proc
{
	// Intersection between the local sub-domain enlarged by the ghost and the contiguous processor
	// sub-domains (External ghost)
	openfpm::vector<::Box<dim,T>> bx;

	// Intersection between the contiguous processor sub-domain enlarged by the ghost with the
	// local sub-domain (Internal ghost)
	openfpm::vector<::Box<dim,T>> nbx;


	// processor
	size_t proc;
};

template<unsigned int dim, typename T>
struct Box_dom
{
	// Intersection between the local sub-domain enlarged by the ghost and the contiguous processor
	// sub-domains (External ghost)
	openfpm::vector_std< Box_sub<dim,T> > ebx;

	// Intersection between the contiguous processor sub-domain enlarged by the ghost with the
	// local sub-domain (Internal ghost)
	openfpm::vector_std< Box_sub<dim,T> > ibx;
};

template<unsigned int dim, typename T>
struct N_box
{
	// id of the processor in the nn_processor list (local processor id)
	size_t id;

	// near processor sub-domains
	typename openfpm::vector<::Box<dim,T>> bx;

	// near processor sector position (or where they live outside the domain)
	openfpm::vector<comb<dim>> pos;

	//! Default constructor
	N_box()
	:id((size_t)-1)
	{};

	//! Copy constructor
	N_box(const N_box<dim,T> & b)
	{
		this->operator=(b);
	}

	//! Copy constructor
	N_box(N_box<dim,T> && b)
	{
		this->operator=(b);
	}

	/*! \brief Copy the element
	 *
	 * \param ele element to copy
	 *
	 */
	N_box<dim,T> & operator=(const N_box<dim,T> & ele)
	{
		id = ele.id;
		bx = ele.bx;
		pos = ele.pos;

		return * this;
	}

	/*! \brief Copy the element
	 *
	 * \param ele element to copy
	 *
	 */
	N_box<dim,T> & operator=(N_box<dim,T> && ele)
	{
		id = ele.id;
		bx.swap(ele.bx);
		pos = ele.pos;

		return * this;
	}

	/*! \brief Compare two N_box object
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

		return bx == ele.bx;
	}

	/*! \brief Compare two N_box object
	 *
	 * \return true if they match
	 *
	 */
	bool operator!=(const N_box<dim,T> & ele) const
	{
		return ! this->operator==(ele);
	}
};

// It store all the boxes of the near processors in a linear array
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
	 */
	bool operator==(const p_box & pb)
	{
		return pb.lc_proc == lc_proc;
	}
};

#endif /* SRC_DECOMPOSITION_COMMON_HPP_ */
