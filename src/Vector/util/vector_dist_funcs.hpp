/*
 * vector_dist_funcs.hpp
 *
 *  Created on: Aug 15, 2018
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_FUNCS_HPP_
#define VECTOR_DIST_FUNCS_HPP_

//! process the particle without properties
struct proc_without_prp
{
	//! process the particle
	template<typename T1, typename T2> inline static void proc(size_t lbl, size_t cnt, size_t id, T1 & v_prp, T2 & m_prp)
	{
		m_prp.get(lbl).set(cnt, v_prp.get(id));
	}
};

//! process the particle without properties
struct proc_without_prp_device
{
	//! process the particle
	template<typename T1, typename T2> __device__ inline static void proc(size_t cnt, size_t id, T1 & v_prp, T2 & m_prp)
	{
		m_prp.set(cnt, v_prp.get(id));
	}
};


/*! \brief Return a valid particle starting from end and tracing back
 *
 * \param end actual opart particle pointer
 * \param end_id actual end particle point
 *
 * \return a valid particle
 *
 */
template<typename m_opart_type>
inline size_t get_end_valid(long int & end, long int & end_id, m_opart_type & m_opart)
{
	end_id--;

	while (end >= 0 && end_id >= 0 && (long int)m_opart.template get<0>(end) == end_id)
	{
		end_id--;
		end--;
	}

	return end_id;
}

//! It process one particle
/*! \brief
 *
 *  \param i particle id
 *  \param end tail of the vector (when a particle is marked to migrate produce an hole, a particle from the tail is used to fill this
 *         hole), and the "real" end go down of one. This parameter contain the real end
 *  \param m_opart for each particle the property 0 contain where the particle has to go, 2 contain to which processor has to go
 *  \param p_map_req it map processor id to request id (consider you have to send to processor 4 7 9) p_map_req contain 4 -> 0,
 *                      7 -> 1, 9 -> 2
 *  \param m_pos sending buffer to fill for position
 *  \param m_prp sending buffer to fill for properties
 *  \param v_pos particle position
 *  \param v_prp particle properties
 *  \param cnt counter for each sending buffer
 *
 */
template<typename proc_class, typename Top,typename Pmr, typename T1, typename T2, typename T3, typename T4>
inline void process_map_particle(size_t i, long int & end, long int & id_end, Top & m_opart, Pmr p_map_req, T1 & m_pos, T2 & m_prp, T3 & v_pos, T4 & v_prp, openfpm::vector<size_t> & cnt)
{
	long int prc_id = m_opart.template get<2>(i);
	size_t id = m_opart.template get<0>(i);

	if (prc_id >= 0)
	{
		size_t lbl = p_map_req.get(prc_id);

		m_pos.get(lbl).set(cnt.get(lbl), v_pos.get(id));
		proc_class::proc(lbl,cnt.get(lbl),id,v_prp,m_prp);

		cnt.get(lbl)++;

		// swap the particle
		// If a particle migrate we have an hole we cover this hole using a particle from the end of the vector
		long int id_valid = get_end_valid(end,id_end,m_opart);

		if (id_valid > 0 && (long int)id < id_valid)
		{
			v_pos.set(id,v_pos.get(id_valid));
			v_prp.set(id,v_prp.get(id_valid));
		}
	}
	else
	{
		// swap the particle
		long int id_valid = get_end_valid(end,id_end,m_opart);

		if (id_valid > 0 && (long int)id < id_valid)
		{
			v_pos.set(id,v_pos.get(id_valid));
			v_prp.set(id,v_prp.get(id_valid));
		}
	}
}


//! It process one particle
template<typename proc_class, typename Top, typename T1, typename T2, typename T3, typename T4>
__device__ inline void process_map_device_particle(unsigned int i, unsigned int offset, Top & m_opart, T1 & m_pos, T2 & m_prp, T3 & v_pos, T4 & v_prp)
{
	size_t id = m_opart.template get<0>(i+offset);

	m_pos.set(i, v_pos.get(id));
	proc_class::proc(i,id,v_prp,m_prp);
}

//! It process one particle
template<typename Top, typename T2, typename T4, unsigned int ... prp>
__device__ inline void process_ghost_device_particle_prp(unsigned int i, unsigned int offset, Top & g_opart, T2 & m_prp, T4 & v_prp)
{
	unsigned int id = g_opart.template get<1>(i+offset) & 0xFFFFFFFF;

	// source object type
	typedef decltype(v_prp.get(id)) encap_src;
	// destination object type
	typedef decltype(m_prp.get(i)) encap_dst;

	// Copy only the selected properties
	object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(id), m_prp.get(i));
}

#endif /* VECTOR_DIST_FUNCS_HPP_ */
