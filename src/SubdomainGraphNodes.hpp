#ifndef SUBDOMAIN_NODES_HPP
#define SUBDOMAIN_NODES_HPP

/* In a decomposition graph each node represent a sub-domain while an edge represent
 * an interaction between sub-domain (it mean that they have to communicate).
 *
 * Here we list the of property that a vertex node can carry with a brief
 * explanation:
 *
 * x = position x of the sub-domain
 * y = position y of the sub-domain
 * z = position z of the sub-domain
 * communication = is the estimated total communication produced by the sub-domain
 * computation = the total computation produced by the sub-domain
 * memory = estimated memory required by the sub-domain
 * id = which processor own this sub-domain
 * sub-id = sub-decomposition where each group of sub-domain is organized in an
 *          hyper-cube
 *
 * Here we list the properties that an edge node can carry with a brief explanation
 *
 * communication = is the estimated communication between sub-domains
 *
 */

/* \brief Sub-domain vertex graph node
 *
 */

struct nm_v
{
	//! The node contain 3 unsigned long integer for communication computation memory and id
	typedef boost::fusion::vector<float[3], size_t, size_t, size_t, size_t, size_t, size_t, size_t> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! type of the positional field
	typedef float s_type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! pos property id in boost::fusion::vector
	static const unsigned int x = 0;
	//! migration property id in boost::fusion::vector
	static const unsigned int migration = 1;
	//! computation property id in boost::fusion::vector
	static const unsigned int computation = 2;
	//! global_id property id in boost::fusion::vector
	static const unsigned int global_id = 3;
	//! id property id in boost::fusion::vector
	static const unsigned int id = 4;
	//! sub_id property id in boost::fusion::vector
	static const unsigned int sub_id = 5;
	//! proc_id property id in boost::fusion::vector
	static const unsigned int proc_id = 6;
	//! fake_v property id in boost::fusion::vector
	static const unsigned int fake_v = 7;

	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 8;

	//! default constructor
	nm_v()
	{

	}

	inline nm_v(const nm_v & p)
	{
		boost::fusion::at_c<0>(data)[0] = boost::fusion::at_c<0>(p.data)[0];
		boost::fusion::at_c<0>(data)[1] = boost::fusion::at_c<0>(p.data)[1];
		boost::fusion::at_c<0>(data)[2] = boost::fusion::at_c<0>(p.data)[2];
		boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
		boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
		boost::fusion::at_c<3>(data) = boost::fusion::at_c<3>(p.data);
		boost::fusion::at_c<4>(data) = boost::fusion::at_c<4>(p.data);
		boost::fusion::at_c<5>(data) = boost::fusion::at_c<5>(p.data);
		boost::fusion::at_c<6>(data) = boost::fusion::at_c<6>(p.data);
		boost::fusion::at_c<7>(data) = boost::fusion::at_c<7>(p.data);
	}

	template<unsigned int dim, typename Mem> inline nm_v(const encapc<dim, nm_v, Mem> & p)
	{
		this->operator=(p);
	}

	template<unsigned int dim, typename Mem> inline nm_v & operator=(const encapc<dim, nm_v, Mem> & p)
	{
		boost::fusion::at_c<0>(data)[0] = p.template get<0>()[0];
		boost::fusion::at_c<0>(data)[1] = p.template get<0>()[1];
		boost::fusion::at_c<0>(data)[2] = p.template get<0>()[2];
		boost::fusion::at_c<1>(data) = p.template get<1>();
		boost::fusion::at_c<2>(data) = p.template get<2>();
		boost::fusion::at_c<3>(data) = p.template get<3>();
		boost::fusion::at_c<4>(data) = p.template get<4>();
		boost::fusion::at_c<5>(data) = p.template get<5>();
		boost::fusion::at_c<6>(data) = p.template get<6>();
		boost::fusion::at_c<7>(data) = p.template get<7>();

		return *this;
	}

	template<unsigned int id> inline auto get() -> decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int id> inline auto get() const -> const decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	static bool noPointers()
	{
		return true;
	}

};

const std::string nm_v::attributes::name[] = { "x", "migration", "computation", "global_id", "id", "sub_id", "proc_id", "fake_v" };

/*! \brief sub-domain edge graph node
 *
 */

struct nm_e
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<size_t, size_t, size_t> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! computation property id in boost::fusion::vector
	static const unsigned int communication = 0;
	static const unsigned int srcgid = 1;
	static const unsigned int dstgid = 2;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 3;

	nm_e()
	{

	}

	template<unsigned int dim, typename Mem> inline nm_e(const encapc<dim, nm_e, Mem> & p)
	{
		boost::fusion::at_c<0>(data) = p.template get<0>();
		boost::fusion::at_c<1>(data) = p.template get<1>();
		boost::fusion::at_c<2>(data) = p.template get<2>();

	}

	template<unsigned int id> inline auto get() -> decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	static bool noPointers()
	{
		return true;
	}
};

const std::string nm_e::attributes::name[] = { "communication", "srcgid", "dstgid" };

/*! \brief Reduced sub-domain vertex graph node
 *
 * It contain only the processor id for each node
 *
 */

struct nm_part_v
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<size_t, size_t> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	typedef float s_type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data

	type data;

	//! partition id in the boost::fusion::vector
	static const unsigned int id = 0;
	//! partition id in the boost::fusion::vector
	static const unsigned int sub_id = 1;

	//! total number of properties
	static const unsigned int max_prop = 2;

	//! default constructor
	nm_part_v()
	{

	}

	template<unsigned int dim, typename Mem> inline nm_part_v(const encapc<dim, nm_part_v, Mem> & p)
	{
		boost::fusion::at_c<0>(data) = p.template get<0>();
		boost::fusion::at_c<1>(data) = p.template get<1>();
	}

};

const std::string nm_part_v::attributes::name[] = { "id", "sub_id" };

/*! \brief Reduced edge graph node
 *
 * It contain only the communication between nodes
 *
 */

struct nm_part_e
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

	//! The data

	type data;

	//! total number of properties
	static const unsigned int max_prop = 0;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};
};

const std::string nm_part_e::attributes::name[] = { "id" };

#endif
