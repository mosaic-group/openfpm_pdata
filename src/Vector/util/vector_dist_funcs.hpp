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

//! It process one particle
template<typename T2, typename T4, unsigned int ... prp>
__device__ inline void process_ghost_device_particle_prp(unsigned int i, unsigned int offset, T2 & m_prp, T4 & v_prp)
{
	unsigned int id = i+offset;

	// source object type
	typedef decltype(v_prp.get(id)) encap_src;
	// destination object type
	typedef decltype(m_prp.get(i)) encap_dst;

	// Copy only the selected properties
	object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(id), m_prp.get(i));
}

template<typename base_type, unsigned int prp>
struct compare_host_device
{
	template<typename St, typename vector_type>
	static bool compare(vector_type & v_prp,St & tol, St & near, bool silent = false)
	{
		bool ret = true;

		// Create a temporal
		openfpm::vector<aggregate<base_type>> tmp;

		tmp.resize(v_prp.size());

		// move host memory to tmp
		auto it = v_prp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			tmp.template get<0>(p) = v_prp.template get<prp>(p);

			++it;
		}

		v_prp.template deviceToHost<prp>();

		// move host memory to tmp
		it = v_prp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			if (fabs(tmp.template get<0>(p) - v_prp.template get<prp>(p)) >= tol && (fabs(tmp.template get<0>(p)) > near && fabs(v_prp.template get<prp>(p)) ) )
			{
				std::cout << "Host and Device buffer differ over set tollerance: " << "Host[" << p << "]="  << tmp.template get<0>(p)
				                                                                   << "  Device[" << p << "]="<< v_prp.template get<0>(p) <<
				                                                                   "  differ more than: " << tol << std::endl;
				ret = false;
			}

			++it;
		}

		//restore
		it = tmp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			v_prp.template get<prp>(p) = tmp.template get<0>(p);

			++it;
		}

		return ret;
	}
};

template<typename base_type,unsigned int N1, unsigned int prp>
struct compare_host_device<Point<N1,base_type>,prp>
{
	template<typename St, typename vector_type>
	static bool compare(vector_type & v_pos,St & tol, St & near, bool silent = false)
	{
		bool ret = true;

		// Create a temporal
		openfpm::vector<Point<N1,base_type>> tmp;

		tmp.resize(v_pos.size());

		// move host memory to tmp
		auto it = v_pos.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			tmp.get(p) = v_pos.get(p);

			++it;
		}

		v_pos.template deviceToHost<prp>();

		// move host memory to tmp
		it = v_pos.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			for (size_t j = 0 ; j < N1 ; j++)
			{
				if (fabs(tmp.template get<0>(p)[j] - v_pos.template get<0>(p)[j]) >= tol && (fabs(tmp.template get<0>(p)[j]) > near && fabs(v_pos.template get<0>(p)[j]) ) )
				{
					std::cout << "Host and Device buffer differ over set tollerance: " << "Host[" << p << "][" << j <<"]="  << tmp.template get<0>(p)[j]
					                                                                 << "  Device[" << p << "][" << j << "]=" << v_pos.template get<0>(p)[j] <<
					                                                                   "  differ more than: " << tol << std::endl;
					ret = false;
				}
			}

			++it;
		}

		//restore
		it = tmp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			v_pos.get(p) = tmp.get(p);

			++it;
		}

		return ret;
	}
};

template<typename base_type,unsigned int N1, unsigned int prp>
struct compare_host_device<base_type[N1],prp>
{
	template<typename St, typename vector_type>
	static bool compare(vector_type & v_prp,St & tol, St & near, bool silent = false)
	{
		bool ret = true;

		// Create a temporal
		openfpm::vector<aggregate<base_type[N1]>> tmp;

		tmp.resize(v_prp.size());

		// move host memory to tmp
		auto it = v_prp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			for (size_t j = 0 ; j < N1 ; j++)
			{
				tmp.template get<0>(p)[j] = v_prp.template get<prp>(p)[j];
			}

			++it;
		}

		v_prp.template deviceToHost<prp>();

		// move host memory to tmp
		it = v_prp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			for (size_t j = 0 ; j < N1 ; j++)
			{
				if (fabs(tmp.template get<0>(p)[j] - v_prp.template get<prp>(p)[j]) >= tol && (fabs(tmp.template get<0>(p)[j]) > near && fabs(v_prp.template get<prp>(p)[j]) ) )
				{
					std::cout << "Host and Device buffer differ over set tollerance: " << "Host[" << p << "]="  << tmp.template get<0>(p)[j]
					                                                                   << "  Device[" << p << "]="<< v_prp.template get<prp>(p)[j] <<
					                                                                   "  differ more than: " << tol << std::endl;
					ret = false;
				}
			}

			++it;
		}

		//restore
		it = v_prp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			for (size_t j = 0 ; j < N1 ; j++)
			{
				v_prp.template get<prp>(p)[j] = tmp.template get<0>(p)[j];
			}

			++it;
		}

		return ret;
	}
};

template<typename base_type,unsigned int N1 , unsigned int N2, unsigned int prp>
struct compare_host_device<base_type[N1][N2],prp>
{
	template<typename St, typename vector_type>
	static bool compare(vector_type & v_prp,St & tol, St & near, bool silent = false)
	{
		bool ret = true;

		// Create a temporal
		openfpm::vector<aggregate<base_type[N1][N2]>> tmp;

		tmp.resize(v_prp.size());

		// move host memory to tmp
		auto it = v_prp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			for (size_t j = 0 ; j < N1 ; j++)
			{
				for (size_t k = 0 ; k < N2 ; k++)
				{
					tmp.template get<0>(p)[j][k] = v_prp.template get<prp>(p)[j][k];
				}
			}

			++it;
		}

		v_prp.template deviceToHost<prp>();

		// move host memory to tmp
		it = v_prp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			for (size_t j = 0 ; j < N1 ; j++)
			{
				for (size_t k = 0 ; k < N2 ; k++)
				{
					if (fabs(tmp.template get<0>(p)[j][k] - v_prp.template get<prp>(p)[j][k]) >= tol && (fabs(tmp.template get<0>(p)[j][k]) > near && fabs(v_prp.template get<prp>(p)[j][k]) ) )
					{
						std::cout << "Host and Device buffer differ over set tollerance: " << "Host[" << p << "]["  << j << "][" << k << "]=" << tmp.template get<0>(p)[j][k]
						                                                                 << "  Device[" << p << "][" << j << "][" << k << "]=" << v_prp.template get<prp>(p)[j][k] << "  differ more than: " << tol << std::endl;

						ret = false;
					}
				}
			}

			++it;
		}

		//restore
		it = v_prp.getIterator();

		while (it.isNext())
		{
			auto p = it.get();

			for (size_t j = 0 ; j < N1 ; j++)
			{
				for (size_t k = 0 ; k < N2 ; k++)
				{
					v_prp.template get<prp>(p)[j][k] = tmp.template get<0>(p)[j][k];
				}
			}

			++it;
		}

		return ret;
	}
};

#endif /* VECTOR_DIST_FUNCS_HPP_ */
