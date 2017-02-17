/*
 * se_class3_vector.hpp
 *
 *  Created on: Feb 11, 2017
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_SE_CLASS3_VECTOR_HPP_
#define SRC_VECTOR_SE_CLASS3_VECTOR_HPP_

#include <iostream>
#include <Space/Shape/Point.hpp>
#include <Vector/map_vector.hpp>
#include <Vector/vector_dist_comm.hpp>
#include <list>

#define SE3_STATUS -2
#define SE3_TYPE -1

enum statuses
{
	CLEAN,
	DIRTY,
	UNINITIALIZED
};

enum sync
{
	SYNC,
	NOTSYNC
};

enum ptype
{
	HALO,
	GHOST,
	INSIDE
};

// Unknown type
template<typename tcheck, bool foundamental>
struct typeCheck
{
	static bool isNan(const tcheck & data)
	{
		return false;
	}

	static bool isInf(const tcheck & data)
	{
		return false;
	}
};

// Unknown type
template<typename tcheck>
struct typeCheck<tcheck,true>
{
	static bool isNan(const tcheck & data)
	{
		return std::isnan(data);
	}

	static bool isInf(const tcheck & data)
	{
		return std::isinf(data);
	}
};

// Array
template<typename tcheck, bool foundamental, unsigned int N1>
struct typeCheck<tcheck[N1], foundamental>
{
	static bool isNan(tcheck (& data)[N1])
	{
		bool nn = false;

		for (size_t i = 0 ; i < N1; i++)
		{
			if (std::isnan(data[i]))
				nn = true;
		}

		return nn;
	}

	static bool isInf(tcheck (& data)[N1])
	{
		bool nn = false;

		for (size_t i = 0 ; i < N1; i++)
		{
			if (std::isinf(data[i]))
				nn = true;
		}

		return nn;
	}
};

// Array2d
template<typename tcheck, bool foundamental, unsigned int N1, unsigned int N2>
struct typeCheck<tcheck[N1][N2], foundamental>
{
	static bool isNan(tcheck (& data)[N1][N2])
	{
		bool nn = false;

		for (size_t i = 0 ; i < N1; i++)
		{
			for (size_t j = 0 ; j < N2; j++)
			{
				if (std::isnan(data[i][j]))
					nn = true;
			}
		}

		return nn;
	}

	static bool isInf(tcheck (& data)[N1][N2])
	{
		bool nn = false;

		for (size_t i = 0 ; i < N1; i++)
		{
			for (size_t j = 0 ; j < N2; j++)
			{
				if (std::isinf(data[i][j]))
					nn = true;
			}
		}

		return nn;
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * property it check that there are not NAN properties
 *
 * \param T boost::fusion::vector
 *
 */
template<typename vector>
struct propCheckNAN
{
	//! Data to check
	vector & data;

	//! Element to check
	size_t id;

	/*! \brief constructor
	 *
	 * \param
	 *
	 */
	inline propCheckNAN(vector & data, size_t id)
	:data(data),id(id)
	{};


	/*!  \brief It call the copy function for each property
	 *
	 * \param t each member
	 *
	 */
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<typename vector::value_type::type,typename boost::mpl::int_<T::value> >::type type_to_check;

		bool snn = typeCheck<type_to_check,std::is_fundamental<type_to_check>::value>::isNan(data.template getProp<T::value>(id));

		if (snn == true)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " error detected NAN in property " << T::value  << std::endl;

			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
	}
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * property it check that there are not NAN properties
 *
 * \param T boost::fusion::vector
 *
 */
template<typename vector>
struct propCheckINF
{
	//! Data to check
	vector & data;

	//! id
	size_t id;


	/*! \brief constructor
	 *
	 * \param
	 *
	 */
	inline propCheckINF(vector & data, size_t id)
	:data(data),id(id)
	{};


	/*!  \brief It call the copy function for each property
	 *
	 * \param t each member
	 *
	 */
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<typename vector::value_type::type,boost::mpl::int_<T::value> >::type type_to_check;

		bool snn = typeCheck<type_to_check,std::is_fundamental<type_to_check>::value>::isInf(data.template getProp<T::value>(id));

		if (snn == true)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " error detected INF in property " << T::value << std::endl;
			ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
		}
	}
};

/*! \brief Return the type of the particle at string level
 *
 * \param type type of the particle
 *
 */
static inline std::string getParticleTypeString(size_t type)
{
	if (type == INSIDE)
		return std::string("INSIDE");
	else if (type == HALO)
		return std::string("HALO");
	else if (type == GHOST)
		return std::string("GHOST");

	return std::string();
}

/*! \brief This class check for inconsistency access
 *
 * \tparam Np number of properties
 *
 */
template<unsigned int Np, unsigned int dim, typename T, typename Decomposition, typename vector>
class se_class3_vector
{
		//! status of the properties
		int sync[2][Np];

		//! number of real properties + POSITION
		static const size_t Np_real = Np+SE3_STATUS+1;

		//! Domain decomposition object
		Decomposition & dec;

		//! Reference to the distributed object
		vector & vd;

		//! temporal buffer
		openfpm::vector<size_t> non_NP;

		bool isLocalHalo(const Point<dim,T> & p)
		{
			for (size_t i = 0; i < dec.getNLocalSub(); i++)
			{
				size_t Nl = dec.getLocalNIGhost(i);

				for (size_t j = 0; j < Nl; j++)
				{
					Box<dim,T> box = dec.getLocalIGhostBox(i, j);

					if (box.isInside(p) == true)
					{
						return true;
					}
				}
			}

			return false;
		}

		/*! \brief Given the position it return the particle type
		 *
		 * \param p position of the particle
		 * \param id of the particle (element id in the vector position)
		 * \param vd reference to the vector
		 *
		 * \return the particle type
		 *
		 */
		size_t getParticleType(const Point<dim,T> & p, const size_t & id, vector & vd)
		{
			size_t type;

			// first we distinguish what is this particle

			if (id > vd.size_local())
				type = GHOST;
			else
			{
				// Use cart decomposition to understand if it is in the halo

				const openfpm::vector<size_t> & vp_id = dec.template ghost_processorID<typename Decomposition::lc_processor_id>(p);

				if (vp_id.size() != 0)
					type = HALO;
				else
				{
					// Check if it is in the HALO inner ghost

					if (isLocalHalo(p) == true)
						type = HALO;
					else
						type = INSIDE;
				}
			}

			return type;
		}

		template<unsigned int ... prp> void create_NNP( const size_t (& gg)[sizeof...(prp)+1] )
		{
			non_NP.clear();

			for (size_t i = 0 ; i < Np_real-1 ; i++)
			{
				bool found = false;
				for (size_t j = 0 ; j < sizeof...(prp) ; j++)
				{
					if (i == gg[j])
					{
						found = true;
						break;
					}
				}

				if (found == false)
					non_NP.add(i);
			}
		}


		std::string getPrpName(size_t i) const
		{
			if (i == Np_real-1)
				return std::string("POSITION");

			return std::to_string(i);
		}

	public:

		//! Constructor all properties are uninitialized
		se_class3_vector(Decomposition & dec, vector & vd)
		:dec(dec),vd(vd)
		{
		}

		template<unsigned int prp> size_t isGhostSync()
		{
			return sync[GHOST][prp];
		}

		void Initialize()
		{
			auto it = vd.getDomainIterator_no_se3();

			while (it.isNext())
			{
				auto p = it.get();

				for (size_t i = 0 ; i < Np_real ; i++)
					vd.template getProp<Np+SE3_STATUS>(p)[i] = UNINITIALIZED;

				vd.template getProp<Np+SE3_TYPE>(p) = INSIDE;

				++it;
			}

			for (size_t i = 0 ; i < Np_real ; i++)
				sync[GHOST][i] = NOTSYNC;
		}

		template <unsigned int ... prp> void ghost_get_pre(size_t opt)
		{
			const size_t gg[sizeof...(prp)+1]  = {prp...};

			create_NNP<prp...>(gg);

			// First check that the ghost are not dirty
			// if they are dirty we are dostroyign information

			auto it = vd.getGhostIterator_no_se3();

			while(it.isNext())
			{
				auto p = it.get();

				for (size_t i = 0 ; i < sizeof...(prp) ; i++)
				{
					if (vd.template getProp<Np+SE3_STATUS>(p)[gg[i]] == DIRTY)
					{
						std::cerr << __FILE__ << ":" << __LINE__ << " Error the ghost has been written and ghost_get will overwrite your changes" << std::endl;
						ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
					}
				}

				if (!(opt & KEEP_PROPERTIES))
				{
					for (size_t i = 0 ; i < non_NP.size() ; i++)
					{
						if (vd.template getProp<Np+SE3_STATUS>(p)[non_NP.get(i)] == DIRTY)
						{
							std::cerr << __FILE__ << ":" << __LINE__ << " Error the it seem that the property=" << getPrpName(non_NP.get(i)) << " has been written and ghost_get will destroy such changes" << std::endl;
							ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
						}
					}
				}

				++it;
			}
		}

		template <unsigned int ... prp> void ghost_get_post(size_t opt)
		{
			const size_t gg[sizeof...(prp)+1]  = {prp...};

			create_NNP<prp...>(gg);

			auto it2 = vd.getGhostIterator_no_se3();

			while(it2.isNext())
			{
				auto p = it2.get();

				for (size_t i = 0 ; i < sizeof...(prp) ; i++)
				{
					if (vd.template getProp<Np+SE3_STATUS>(p)[gg[i]] == DIRTY)
						vd.template getProp<Np+SE3_STATUS>(p)[gg[i]] = CLEAN;
				}

				vd.template getProp<Np+SE3_TYPE>(p) = GHOST;

				++it2;
			}

			if (!(opt & KEEP_PROPERTIES))
			{
				for (size_t i = 0 ; i < non_NP.size() ; i++)
					sync[GHOST][gg[i]] = NOTSYNC;

				auto it = vd.getGhostIterator_no_se3();

				while (it.isNext() == true)
				{
					auto p = it.get();

					for (size_t i = 0 ; i < non_NP.size() ; i++)
						vd.template getProp<Np+SE3_STATUS>(p)[non_NP.get(i)] = UNINITIALIZED;

					++it;
				}
			}

			// We notify that the ghost are in sync

			for (size_t i = 0 ; i < sizeof...(prp) ; i++)
				sync[GHOST][gg[i]] = SYNC;

			if (!(opt & NO_POSITION))
				sync[GHOST][Np_real] = SYNC;

			if (!(opt & KEEP_PROPERTIES))
			{
				for (size_t i = 0 ; i < non_NP.size() ; i++)
					sync[GHOST][non_NP.get(i)] = NOTSYNC;
			}
		}

		template <unsigned int ... prp> void ghost_put()
		{
			const size_t gg[sizeof...(prp)]  = {prp...};

			auto it = vd.getDomainIterator_no_se3();

			while(it.isNext())
			{
				auto p = it.get();

				if (vd.template getProp<Np+SE3_TYPE>(p) == INSIDE)
				{
					++it;
					continue;
				}

				for (size_t i = 0 ; i < sizeof...(prp) ; i++)
				{
					if (vd.template getProp<Np+SE3_STATUS>(p)[gg[i]] == UNINITIALIZED)
					{
						std::cerr << __FILE__ << ":" << __LINE__ << " error it seem that you are sending at least in part uninitialized data with ghost put " << std::endl;
						ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
					}
				}

				++it;
			}

			for (size_t i = 0 ; i < sizeof...(prp) ; i++)
			{
				sync[HALO][gg[i]] = SYNC;
				sync[HALO][gg[i]] = CLEAN;
			}

			// Ghost has been merged make ghost clean
			auto it2 = vd.getGhostIterator_no_se3();

			while(it2.isNext())
			{
				auto p = it2.get();

				for (size_t i = 0 ; i < sizeof...(prp) ; i++)
					vd.template getProp<Np+SE3_STATUS>(p)[gg[i]] = CLEAN;

				++it2;
			}
		}

		void getIterator() const
		{
			auto it = vd.getDomainIterator_no_se3();

			while(it.isNext())
			{
				auto p = it.get();

				for (size_t j = 0 ; j < Np ; j++)
				{
					if (vd.template getProp<Np+SE3_STATUS>(p)[j] == DIRTY)
						vd.template getProp<Np+SE3_STATUS>(p)[j] = CLEAN;
				}

#ifdef CHECK_FOR_POSINF

				if ( std::isinf(vd.getPos(p)[0]) || std::isinf(vd.getPos(p)[1]) || std::isinf(vd.getPos(p)[2]) )
				{
					std::cerr << __FILE__ << ":" << __LINE__ << " error detected INF in position for particle p=" << p.getKey() << " of type=" << getParticleTypeString(vd.template getProp<Np+SE3_TYPE>(p)) << std::endl;
					ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
				}

#endif

#ifdef CHECKFOR_POSNAN

				if ( std::isnan(vd.getPos(p)[0]) || std::isnan(vd.getPos(p)[1]) || std::isnan(vd.getPos(p)[2]) )
				{
					std::cerr << __FILE__ << ":" << __LINE__ << " error detected NAN in position for particle p=" << p.getKey() << " of type=" << getParticleTypeString(vd.template getProp<Np+SE3_TYPE>(p)) << std::endl;
					ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
				}

#endif

#ifdef CHECKFOR_PROPINF

				{
					propCheckINF<vector> checker(vd,p.getKey());

					boost::mpl::for_each_ref< boost::mpl::range_c<int,0, Np_real-1 > > (checker);
				}

#endif

#ifdef CHECKFOR_PROPNAN

				{
					propCheckNAN<vector> checker(vd,p.getKey());

					boost::mpl::for_each_ref< boost::mpl::range_c<int,0, Np_real-1 > >(checker);
				}

#endif

				++it;
			}
		}

		void map_pre()
		{
			auto it = vd.getGhostIterator_no_se3();

			while (it.isNext() == true)
			{
				auto p = it.get();

				for (size_t j = 0 ; j < Np ; j++)
				{
					if (vd.template getProp<Np+SE3_STATUS>(p)[j] == DIRTY)
					{
						std::cerr << __FILE__ << ":" << __LINE__ << " error it seem that ghost has been filled with information that we are going to destroy with the map call " << std::endl;
					}
				}

				++it;
			}
		}

		void map_post()
		{
			for (size_t j = 0 ; j < Np ; j++)
			{

				sync[GHOST][j] = NOTSYNC;
			}

			auto it = vd.getDomainIterator_no_se3();

			while (it.isNext() == true)
			{
				auto p = it.get();

				for (size_t j = 0 ; j < Np ; j++)
				{
					Point<vector::dims,typename vector::stype> xp = vd.getPos(p);

					vd.template getProp<Np+SE3_TYPE>(p) = getParticleType(xp,p.getKey(),vd);
				}

				++it;
			}
		}

		template<unsigned int prp> void read(const vector & vd, size_t p) const
		{
			if (vd.template getProp<Np+SE3_STATUS>(p)[prp] == UNINITIALIZED)
			{
				std::stringstream str;
				std::string type_str = getParticleTypeString(vd.template getProp<Np+SE3_TYPE>(p));

				if (prp == Np_real)
					str << __FILE__ << ":" << __LINE__ << " Error you are reading the particle " << p << " of type " << type_str << " the position. But it result to be uninitialized" << std::endl;
				else
					str << __FILE__ << ":" << __LINE__ << " Error you are reading from the particle " << p << " of type " << type_str << " the property=" << getPrpName(prp) << ". But it result to be uninitialized" << std::endl;

				// It is an error read from an uninitialized property
				std::cerr << str.str();
				ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
			}

			if (vd.template getProp<Np+SE3_STATUS>(p)[prp] == DIRTY)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Warning you are reading from an particle that has been changed already in the same cycle" << std::endl;
				ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
			}

			if (vd.template getProp<Np+SE3_TYPE>(p) == GHOST || vd.template getProp<Np+SE3_TYPE>(p) == HALO)
			{
				// if we read from the ghost we have to ensure that the ghost is in
				// sync in particular that the state of the halo is CLEAR

				if (sync[vd.template getProp<Np+SE3_TYPE>(p)][prp] != SYNC)
				{
					std::cerr << __FILE__ << ":" << __LINE__ << " Error it seem that you are reading from a ghost the property=" << getPrpName(prp) << " but it seem it is changed from the last ghost_get. It seems that is missing a ghost_get" << std::endl;
					ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
				}
			}
		}

		template<unsigned int prp> void write(vector & vd, size_t p)
		{
			vd.template getProp<Np+SE3_STATUS>(p)[prp] = DIRTY;
			if (p >= vd.size_local())
				vd.get_se_class3().template setHaloOutSync<prp>();
			else
			{
				if (vd.template getProp<Np+SE3_TYPE>(p) == HALO)
					vd.get_se_class3().template setGhostOutSync<prp>();
			}

		}

		//! Copy operator
		se_class3_vector<Np,dim,T,Decomposition,vector> & operator=(const se_class3_vector<Np,dim,T,Decomposition,vector> & se3)
		{
			for (size_t i = 0 ; i < Np ; i++)
			{
				sync[0][i] = se3.sync[0][i];
				sync[1][i] = se3.sync[1][i];
			}

			dec = se3.dec;
			vd = se3.vd;

			return *this;
		}

		template<unsigned int prp> void setHaloOutSync()
		{
			for (size_t i = 0 ; i < Np ; i++)
				sync[HALO][prp] = NOTSYNC;
		}

		template<unsigned int prp> void setGhostOutSync()
		{
			for (size_t i = 0 ; i < Np ; i++)
				sync[GHOST][prp] = NOTSYNC;
		}

		void getNN()
		{
			if (sync[GHOST][Np_real] == NOTSYNC)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " Error you are trying to get a Cell-list or Verlet-list without having the ghost synchronized in position please use ghost_get before" << std::endl;
				ACTION_ON_ERROR(VECTOR_DIST_ERROR_OBJECT);
			}
		}
};


#endif /* SRC_VECTOR_SE_CLASS3_VECTOR_HPP_ */
