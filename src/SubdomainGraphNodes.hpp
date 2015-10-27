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
	typedef boost::fusion::vector<float,float,float,size_t,size_t,size_t,size_t,long int> type;

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

	//! computation property id in boost::fusion::vector
	static const unsigned int x = 0;
	//! computation property id in boost::fusion::vector
	static const unsigned int y = 1;
	//! memory property id in boost::fusion::vector
	static const unsigned int z = 2;
	//! computation property id in boost::fusion::vector
	static const unsigned int communication = 3;
	//! computation property id in boost::fusion::vector
	static const unsigned int computation = 4;
	//! memory property id in boost::fusion::vector
	static const unsigned int memory = 5;
	//! memory property id in boost::fusion::vector
	static const unsigned int id = 6;
	//! memory property sub_id in boost::fusion::vector
	static const unsigned int sub_id = 7;

	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 8;
    
    //! default constructor
    nm_v(){
        
    }
    
    template <unsigned int dim, typename Mem> inline nm_v(const encapc<dim,nm_v,Mem> & p)
    {
        boost::fusion::at_c<0>(data) = p.template get<0>();
        boost::fusion::at_c<1>(data) = p.template get<1>();
        boost::fusion::at_c<2>(data) = p.template get<2>();
        boost::fusion::at_c<3>(data) = p.template get<3>();
        boost::fusion::at_c<4>(data) = p.template get<4>();
        boost::fusion::at_c<5>(data) = p.template get<5>();
        boost::fusion::at_c<6>(data) = p.template get<6>();
        boost::fusion::at_c<7>(data) = p.template get<7>();
    }
    
    
};

const std::string nm_v::attributes::name[] = {"x","y","z","communication","computation","memory","id","sub_id"};

/*! \brief sub-domain edge graph node
 *
 */

struct nm_e
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<size_t> type;

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
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 1;
};

const std::string nm_e::attributes::name[] = {"communication"};

/*! \brief Reduced sub-domain vertex graph node
 *
 * It contain only the processor id for each node
 *
 */

struct nm_part_v
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<size_t,size_t> type;

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
    nm_part_v(){
        
    }
    
    template <unsigned int dim, typename Mem> inline nm_part_v(const encapc<dim,nm_part_v,Mem> & p)
    {
        boost::fusion::at_c<0>(data) = p.template get<0>();
        boost::fusion::at_c<1>(data) = p.template get<1>();
    }
    
};



const std::string nm_part_v::attributes::name[] = {"id","sub_id"};

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

const std::string nm_part_e::attributes::name[] = {"id"};

#endif
