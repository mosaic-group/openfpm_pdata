/*
 * dist_map_graph.hpp
 *
 *  Created on: Dec 10, 2015
 *      Author: Antonio Leo
 *
 *
 *  The distributed graph, is a graph distributed across processors, each processor store part of the graph
 *
 * ## Dictionary
 *
 * * local vertex id, is the index of the vertex in the local graph
 * * vertex id is the "unique" index in the distribution vector (vtxdist)
 * * global id is the "unique" topological id of the vertex, it never change
 *
 *
 *
 *  Graph structure that store a CSR graph format
 *
 *  This Graph format is suppose to have a list of vertex that store an index x that indicate where
 *  in the adjacency list we have the list of all the neighborhood vertex, plus an adjacency list
 *  that store all the neighborhood for each vertex:
 *
 *  In reality inside DistGraph_CSR
 *
 *  VertexList store the Vertex index that indicate the end of the neighborhood list,
 *  the start is indicated by i * v_slot (id of the vertex, v_slot maximum number of
 *  neighborhood for a vertex)
 *
 *  EdgeList store for each vertex at position i * v_slot (id of the vertex, v_slot maximum
 *  number of neighborhood for a vertex) the list of all the neighborhood of the vertex i
 *
 *  Example
 *
 *  Suppose an undirected graph of 4 vertex
 *
 *  1 -> 2 3 4
 *  2 -> 1
 *  3 -> 4 1
 *  4 -> 1 3
 *
 *  suppose that v_slot is 4 ( each vertex has less tha 4 neighborhood  )
 *
 *  we will have
 *
 *  Vertex list 3 5 10 14
 *  Edge list   2 3 4 0 1 0 0 0 4 1 0 0 1 3 0 0
 *
 *  Vertex properties and edge properties are stored in a separate structure
 *
 * ## Dictionary
 *
 * The distributed
 *
 * * local vertex id is the index of the vertex in
 *
 */

#ifndef DIST_MAP_GRAPH_HPP_
#define DIST_MAP_GRAPH_HPP_

#include "Vector/map_vector.hpp"
#include "Graph/map_graph.hpp"
#include <unordered_map>
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"
#include "VCluster/VCluster.hpp"

#define NO_EDGE -1
#define DIST_GRAPH_ERROR 7001

template<typename V, typename E,
		 typename Memory,
		 typename layout_v,
		 typename layout_e,
		 template<typename> class layout_v_base,
		 template<typename> class layout_e_base,
		 typename grow_p>
class DistGraph_CSR;

class v_info
{
public:
	typedef boost::fusion::vector<size_t, size_t> type;

	type data;

	static const unsigned int id = 0;
	static const unsigned int gid = 1;
	static const unsigned int max_prop = 2;

	v_info()
	{
	}

	inline void setid(size_t id_)
	{
		boost::fusion::at_c<0>(data) = id_;
	}

	inline void setgid(size_t gid_)
	{
		boost::fusion::at_c<1>(data) = gid_;
	}

	template<unsigned int id> inline auto get() -> decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int id> inline auto get() const -> const decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int dim, typename Mem> inline v_info(const encapc<dim, v_info, Mem> & p)
	{
		this->operator=(p);
	}

	template<unsigned int dim, typename Mem> inline v_info & operator=(const encapc<dim, v_info, Mem> & p)
	{
		boost::fusion::at_c<0>(data) = p.template get<0>();
		boost::fusion::at_c<1>(data) = p.template get<1>();

		return *this;
	}

	static bool noPointers()
	{
		return true;
	}
};

class e_info
{
public:
	typedef boost::fusion::vector<size_t, size_t> type;

	type data;

	static const unsigned int sgid = 0;
	static const unsigned int dgid = 1;
	static const unsigned int max_prop = 2;

	e_info()
	{
	}

	template<unsigned int id> inline auto get() -> decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int id> inline auto get() const -> const decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int dim, typename Mem> inline e_info(const encapc<dim, e_info, Mem> & p)
	{
		this->operator=(p);
	}

	template<unsigned int dim, typename Mem> inline e_info & operator=(const encapc<dim, e_info, Mem> & p)
	{
		boost::fusion::at_c<0>(data) = p.template get<0>();
		boost::fusion::at_c<1>(data) = p.template get<1>();

		return *this;
	}

	static bool noPointers()
	{
		return true;
	}
};

/*! \brief Structure that store a graph in CSR format or basically in compressed adjacency matrix format
 *
 * \param V each vertex will encapsulate have this type
 * \param E each edge will encapsulate this type
 * \param device Type of device / basicaly it select the layout
 *        for device_cpu is (x_1, p1_1, p2_1, p3_1 ....), ... ( x_n, p1_1, p2_1, p3_1, ...)
 *        for device_gpu is (x_1, ... , x_n) ... (p1_n, ... pn_n)
 *        where x_1 is the index where it end the list of the neighborhood list and pj_k is
 *        the property j for the vertex j. Basically in the first case one array will store
 *        index and property of each vertex, in the second case several array will store
 *        index and property
 *
 * \param VertexList structure that store the list of Vertex
 * \param EdgeList structure that store the list of edge
 *
 * \warning This graph is suitable only when we know the graph structure and we build
 * the graph adding vertexes and edges, removing vertex and edge is EXTREMLY expensive
 *
 * ### Example on creating graph, moving vertices and redistribution
 * \snippet dist_graph_unit_tests.hpp
 *
 */

template<typename V, typename E = no_edge,
		 typename Memory = HeapMemory,
		 typename layout_v = typename memory_traits_lin<V>::type,
		 typename layout_e = typename memory_traits_lin<E>::type,
		 template <typename> class layout_v_base = memory_traits_lin,
		 template <typename> class layout_e_base = memory_traits_lin,
		typename grow_p = openfpm::grow_policy_double>
class DistGraph_CSR
{
	//! Vcluster communication object
	Vcluster & vcl;

	//! Distribution vector
	openfpm::vector<idx_t> vtxdist;

	//! Fixed distribution vector, it never changes, it maintains always the first decomposition and topology
	openfpm::vector<idx_t> fvtxdist;

	//! number of slot per vertex
	size_t v_slot;

	//! Structure that store the vertex properties
	openfpm::vector<V, Memory, layout_v,layout_v_base,grow_p, openfpm::vect_isel<V>::value> v;

	//! Structure that store the vertex id and global id
	openfpm::vector<v_info, Memory, typename memory_traits_lin<v_info>::type, memory_traits_lin, grow_p, openfpm::vect_isel<v_info>::value> v_m;

	//! Structure that store the number of adjacent vertex in e_l for each vertex
	openfpm::vector<size_t, Memory, typename layout_v_base<size_t>::type, layout_v_base, grow_p, openfpm::vect_isel<size_t>::value> v_l;

	//! Structure that store the edge properties
	openfpm::vector<E, Memory, layout_e, layout_e_base, grow_p, openfpm::vect_isel<E>::value> e;

	//! Structure that store the edge properties
	openfpm::vector<e_info, Memory, typename layout_e_base<e_info>::type, layout_e_base, grow_p, openfpm::vect_isel<e_info>::value> e_m;

	//! Structure that store for each vertex the adjacent the vertex id and edge id (for property into e)
	openfpm::vector<e_map, Memory, typename memory_traits_lin<e_map>::type, layout_e_base, grow_p, openfpm::vect_isel<e_map>::value> e_l;

	//! invalid edge element, when a function try to create an in valid edge this object is returned
	openfpm::vector<E, Memory, layout_e, layout_e_base, grow_p, openfpm::vect_isel<E>::value> e_invalid;

	//! Map to access to the global vertex id given the vertex id
	std::unordered_map<size_t, size_t> id2glb;

	//! Map to access the vertex id given the global vertex id
	std::unordered_map<size_t, size_t> glb2id;

	//! Map to access the local vertex id given the global one
	std::unordered_map<size_t, size_t> glb2loc;

	//! Struct containing the (sub)graph to send
	typedef struct
	{
		//! vertex send buffer
		openfpm::vector<V> send_v;
		//! vertex info send buffer
		openfpm::vector<v_info> send_v_m;
		//! edge send buffer
		openfpm::vector<E> send_e;
		//! edge info send buffer
		openfpm::vector<e_info> send_e_m;
		//! For each edge contain the child vertex id
		openfpm::vector<size_t> send_el;
		//! For each vertex contain the number of children
		openfpm::vector<size_t> send_es;
		//! Indicates if the pack is empty or not
		bool isEmpty = true;
	} SendGraphPack;

	//! Pack storing that data to send to other processors
	openfpm::vector<SendGraphPack> sgp;

	//! Array containing the sent vertices and that will be deleted from the graph
	openfpm::vector<size_t> v_td;

	//! Structure needed to get vertex position by global id
	typedef struct
	{
		//! vertex id
		size_t id;
		// processor containing the vertex
		size_t pid;
	} GlobalVInfo;

	//!TODO update description from pdf
	// Map of GlobalVInfo containing informations of vertices of the INITIAL distribution contained in this processor
	// ex. if this will contain the first 4 vertices of the distribution (0,1,2,3) it will maintain informations only about these vertices
	// The key is the vertex global id
	std::unordered_map<size_t, GlobalVInfo> glbi_map;

	//! Queue of vertex requests
	openfpm::vector<openfpm::vector<size_t>> vr_queue;

	//! Map containing the ghost vertices of this graph, if bool is false the ghost will be deleted in the next vertices exchange
	std::unordered_map<size_t, bool> ghs_map;

	//! Structure to store a add request of an edge
	typedef struct
	{
		//! source vertex
		size_t v1;
		//! target vertex
		size_t v2;
		//! source vertex global index
		size_t v1n;
		//! destination vertex global index
		size_t v2n;
	} EdgeReq;

	//! Queue of requests to add edges
	openfpm::vector<EdgeReq> e_queue;

	/*! \brief add edge on the graph
	 *
	 * add edge on the graph
	 *
	 * \param v1 start vertex
	 * \param v2 end vertex
	 *
	 * \return the index of the edge created
	 *
	 */
	template<typename CheckPolicy = NoCheck> inline size_t addEdge_(size_t v1, size_t v2)
	{
		// If v1 and v2 does not satisfy some criteria return
		if (CheckPolicy::valid(v1, v.size()) == false)
			return NO_EDGE;
		if (CheckPolicy::valid(v2, v.size()) == false)
			return NO_EDGE;

		// get the number of adjacent vertex
		size_t id_x_end = v_l.template get<0>(v1);

#ifdef DEBUG

		// Check that the edge does not already exist

		for (size_t s = 0; s < id_x_end; s++)
		{
			if (e_l.template get<e_map::vid>(v1 * v_slot + s) == v2)
			{
				std::cerr << "Error graph: the edge already exist\n";
			}
		}
#endif

		// Check if there is space for another edge

		if (id_x_end >= v_slot)
		{
			// Unfortunately there is not space we need to reallocate memory
			// Reallocate with double slot

			// Create an new Graph
			DistGraph_CSR<V, E> g_new(2 * v_slot, v.size());

			// Copy the graph
			for (size_t i = 0; i < v.size(); i++)
			{
				// copy the property from the old graph
				g_new.v.set(i, v, 2 * i);
			}

			// swap the new graph with the old one
			swap(g_new);
		}

		// Here we are sure than v and e has enough slots to store a new edge
		// Check that e_l has enough space to store new edge
		// should be always e.size == e_l.size

		if (id_x_end >= e_l.size())
		{
			// Resize the basic structure
			e_l.resize(v.size() * v_slot);
		}

		// add in e_l the adjacent vertex for v1 and fill the edge id
		e_l.template get<e_map::vid>(v1 * v_slot + id_x_end) = v2;
		e_l.template get<e_map::eid>(v1 * v_slot + id_x_end) = e.size();

		// add an empty edge
		e.resize(e.size() + 1);
		e_m.resize(e_m.size() + 1);

		// Increment the ending point
		++v_l.template get<0>(v1);

		// return the created edge
		return e.size() - 1;
	}

	/*! \brief Callback to set the size of the receiving vector
	 *
	 * \param msg_i Index of the message
	 * \param total_msg Total number of messages
	 * \param total_p Total number of processors to communicate with
	 * \param i Processor id
	 * \param ri Request id
	 * \param ptr Void pointer parameter for additional data to pass to the call-back
	 */
	static void * gr_receive(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, size_t tag, void * ptr)
	{
		openfpm::vector<HeapMemory> *v = static_cast<openfpm::vector<HeapMemory> *>(ptr);

		v->get(i).allocate(msg_i);

		return v->get(i).getPointer();

	}

	/*! \brief Callback of the sendrecv to set the size of the array received
	 *
	 * \param msg_i Index of the message
	 * \param total_msg Total number of messages
	 * \param total_p Total number of processors to communicate with
	 * \param i Processor id
	 * \param ri Request id
	 * \param ptr Void pointer parameter for additional data to pass to the call-back
	 */
	static void * on_receive(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, size_t tag, void * ptr)
	{
		openfpm::vector<openfpm::vector<size_t>> *v = static_cast<openfpm::vector<openfpm::vector<size_t>> *>(ptr);

		v->get(i).resize(msg_i / sizeof(size_t));

		return &(v->get(i).get(0));
	}

	/*! \brief Init communication structures
	 *
	 */
	void resetExchange()
	{
		if (sgp.size() == vcl.getProcessingUnits())
		{
			for (size_t p = 0; p < vcl.getProcessingUnits(); ++p)
			{
				sgp.get(p).send_v.clear();
				sgp.get(p).send_v_m.clear();
				sgp.get(p).send_e.clear();
				sgp.get(p).send_e_m.clear();
				sgp.get(p).send_es.clear();
				sgp.get(p).send_el.clear();
				sgp.get(p).isEmpty = true;
			}
		}
		else
		{
			sgp.resize(vcl.getProcessingUnits());

			for (size_t p = 0; p < vcl.getProcessingUnits(); ++p)
			{
				openfpm::vector<V> s_v;
				openfpm::vector<v_info> s_v_m;
				openfpm::vector<E> s_e;
				openfpm::vector<e_info> s_e_m;
				openfpm::vector<size_t> s_el;
				openfpm::vector<size_t> s_es;

				sgp.get(p).send_v = s_v;
				sgp.get(p).send_v_m = s_v_m;
				sgp.get(p).send_e = s_e;
				sgp.get(p).send_e_m = s_e_m;
				sgp.get(p).send_es = s_es;
				sgp.get(p).send_el = s_el;
				sgp.get(p).isEmpty = true;
			}
		}
	}

	/*! \brief Remove from this graph the vertices that have been sent
	 *
	 */
	void deleteMovedVertices()
	{
		// Prepare new graph
		DistGraph_CSR recv_g;

		// Reset maps
		glb2loc.clear();
		glb2id.clear();
		id2glb.clear();

		size_t local_i = 0;

		// Add previous vertices without the sent ones
		for (size_t i = 0; i < this->getNVertex(); ++i)
		{
			if (!isToDelete(i))
			{
				recv_g.add_vertex(this->vertex(i), this->getVertexId(i), this->getVertexGlobalId(i));

				for (size_t j = 0; j < this->getNChilds(i); j++)
				{
					recv_g.addEdge(local_i, this->getChild(i, j), this->getChildEdge(i, j), this->getChildInfo(i, j));
				}
				++local_i;
			}
		}

		recv_g.fvtxdist = fvtxdist;

		// Swap temporary graph with the main one
		swap(recv_g);

		// Clear vertex to delete array
		v_td.clear();
	}

	/*! \brief Check it the vertex i must be deleted or not
	 *
	 * \param i vertex
	 * \return true if to delete, false otherwise
	 */
	bool isToDelete(size_t i)
	{

		for (size_t j = 0; j < v_td.size(); ++j)
		{
			if (i == v_td.get(j))
				return true;
		}
		return false;
	}

	/*! \brief Get the processor of the the given vertex id, CAN be used BEFORE re-mapping starts
	 *
	 * \param v id of the vertex
	 * \return process id
	 */
	size_t getVProcessor(size_t v)
	{
		for (size_t i = 1; i < vtxdist.size() - 1; ++i)
		{
			if (v < vtxdist.get(i))
			{
				return i - 1;
			}
		}
		return vcl.getProcessingUnits() - 1;
	}

	/*! \brief Send and receive vertices and update current graph
	 *
	 * \tparam Remove the sent sub-graph
	 *
	 */
	template<bool addAsGhosts>
	void exchangeVertices()
	{
		// If the exchange is not to retrieve ghost vertices delete the vertices this processor is sending
		if (!addAsGhosts)
			deleteMovedVertices();

		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> size;
		openfpm::vector<void *> ptr;
		openfpm::vector<HeapMemory> packs(vcl.getProcessingUnits());

		// Total number of vertex to send
		size_t nvts = 0;
		for (size_t i = 0; i < vcl.getProcessingUnits(); i++)
		{
			nvts += sgp.get(i).send_v.size();
		}

		// For each processor
		for (size_t i = 0; i < vcl.getProcessingUnits(); i++)
		{
			// if nothing to send continue
			if (sgp.get(i).isEmpty)
				continue;

			// process to communicate with TODO remove pc
			size_t pc = i;

			size_t vp_size = sgp.get(pc).send_v.size();

			size_t req = 0;

			// prepare slot for number of vertices
			Packer<size_t, HeapMemory>::packRequest(req);

			for (size_t j = 0; j < vp_size; j++)
			{
				// prepare slot for vertex
				Packer<V, HeapMemory>::packRequest(req);

				// prepare slot info for vertex
				Packer<v_info, HeapMemory>::packRequest(req);

				// prepare slot for the number of children
				Packer<size_t, HeapMemory>::packRequest(req);

				// prepare slots for the children
				for (size_t k = 0; k < sgp.get(pc).send_es.get(j); k++)
				{
					// prepare slot for edge
					Packer<E, HeapMemory>::packRequest(req);

					// prepare slot for edge info
					Packer<e_info, HeapMemory>::packRequest(req);

					// prepare slot for edge target id
					Packer<size_t, HeapMemory>::packRequest(req);
				}
			}

			// allocate the memory
			HeapMemory & pmem = *(new HeapMemory());
//			pmem.allocate(req);
			ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(req, pmem));

			mem.incRef();

			Pack_stat sts;
			size_t e_it = 0;

			// Pack total size
			Packer<size_t, HeapMemory>::pack(mem, vp_size, sts);

			for (size_t j = 0; j < vp_size; j++)
			{
				// Pack the vertex
				Packer<decltype(sgp.get(pc).send_v.get(0)), HeapMemory>::pack(mem, sgp.get(pc).send_v.get(j), sts);

				// Pack the vertex info
				Packer<decltype(sgp.get(pc).send_v_m.get(0)), HeapMemory>::pack(mem, sgp.get(pc).send_v_m.get(j), sts);

				// Pack size of the children
				Packer<size_t, HeapMemory>::pack(mem, sgp.get(pc).send_es.get(j), sts);

				// Pack children
				for (size_t k = 0; k < sgp.get(pc).send_es.get(j); k++)
				{
					// Pack the edge
					Packer<decltype(sgp.get(pc).send_e.get(0)), HeapMemory>::pack(mem, sgp.get(pc).send_e.get(e_it), sts);

					// Pack the edge info
					Packer<decltype(sgp.get(pc).send_e_m.get(0)), HeapMemory>::pack(mem, sgp.get(pc).send_e_m.get(e_it), sts);

					// Pack the edge target id
					Packer<size_t, HeapMemory>::pack(mem, sgp.get(pc).send_el.get(e_it), sts);

					++e_it;
				}
			}

			prc.add(i);
			size.add(pmem.size());
			ptr.add(mem.getPointerBase());
		}

		// Exchange informations through processors
		vcl.sendrecvMultipleMessagesNBX(prc.size(), (size_t *)size.getPointer(), (size_t *)prc.getPointer(), (void **)ptr.getPointer(), gr_receive, &packs, NONE);

		for (size_t i = 0; i < vcl.getProcessingUnits(); i++)
		{

			if (i != vcl.getProcessUnitID() && packs.get(i).size() > 0)
			{
				Unpack_stat ps;

				ExtPreAlloc<HeapMemory> mem(packs.get(i).size(), packs.get(i));

				// unpack total number of vertex
				size_t r_size;
				Unpacker<size_t, HeapMemory>::unpack(mem, r_size, ps);

				// take previous last item
				size_t prev = getNVertex();

				for (size_t j = prev; j < prev + r_size; j++)
				{
					// unpack the vertex
					V v_n;
					Unpacker<V, HeapMemory>::unpack(mem, v_n, ps);

					v_info vm;
					Unpacker<v_info, HeapMemory>::unpack(mem, vm, ps);

					add_vertex(v_n, vm.template get<v_info::id>(), vm.template get<v_info::gid>());

					if (addAsGhosts)
						ghs_map.insert( { vm.template get<v_info::gid>(), false });

					// unpack size of children
					size_t s;
					Unpacker<size_t, HeapMemory>::unpack(mem, s, ps);

					// prepare slots for the children
					for (size_t k = 0; k < s; k++)
					{
						// unpack edge
						E e_n;
						Unpacker<E, HeapMemory>::unpack(mem, e_n, ps);

						// unpack edge
						e_info e_i;
						Unpacker<e_info, HeapMemory>::unpack(mem, e_i, ps);

						// unpack vertex id of the edge target
						size_t el_n;
						Unpacker<size_t, HeapMemory>::unpack(mem, el_n, ps);

						// add the edge //HERE ERROR modify to add globals
						addEdge(j, el_n, e_n, e_i);
					}
				}
			}
		}

		// After the exchange reset all the structures needed for it
		resetExchange();
	}

	/*! \brief Update the distribution vector vtxdist
	 *
	 */
	void updateVtxdist()
	{
		// Prepare vtxdist
		vtxdist.resize(vcl.getProcessingUnits() + 1);

		// Set first element to 0, always the same
		vtxdist.get(0) = 0;

		// Prepare array that will contain the sizes of all the graphs
		openfpm::vector<size_t> recv(vcl.getProcessingUnits());

		// Set the local size
		size_t tv = getNVertex();

		// Sent and receive the size of each subgraph
		vcl.allGather(tv, recv);
		vcl.execute();

		// Set the value for this processor
		recv.get(vcl.getProcessUnitID()) = getNVertex();

		// Update vtxdist
		for (size_t i = 1; i <= recv.size(); ++i)
		{
			vtxdist.get(i) = recv.get(i - 1) + vtxdist.get(i - 1);
		}
	}

	/*! \brief Re-map received vertices in order to have ordered vertex ids
	 *
	 */
	void remap()
	{
		size_t p_id = vcl.getProcessUnitID();

		typedef struct
		{
			// new vertex id
			size_t id;
			// processor rank that contain the vertex
			size_t pid;
		} IdnProc;

		// Map that will contain the couples to update the global info map in this processor
		// The key is the (old vertex id)
		std::unordered_map<size_t, IdnProc> on_toup;

		// For each processor old, new couples
		openfpm::vector<openfpm::vector<size_t>> on_info(vcl.getProcessingUnits());

		std::map<size_t, size_t> old_glob2loc(glb2loc.begin(), glb2loc.end());
		size_t j = vtxdist.get(p_id);
		size_t i = 0, k = 0;

		// Reset maps of vertices ids
		id2glb.clear();
		glb2id.clear();
		glb2loc.clear();

		// Fix sending couples gid, newid and remove on_toup updating glbi_map here and after receive
		for (auto it : old_glob2loc)
		{
			// The couple to the final map, needed to update the vertices in this sub-graph
			IdnProc nidpid = { j, p_id };
			on_toup.insert( { v_m.get(it.second).template get<v_info::id>(), nidpid });

			// fill the re-mapping information for each processors that need it
			on_info.get(getInfoProc(it.first)).add(v_m.get(it.second).template get<v_info::id>());
			on_info.get(getInfoProc(it.first)).add(j);

			// Re-map the vertex
			map_v(j, v_m.get(it.second).template get<v_info::gid>(), it.second);

			j++;
			i++;
			k += 2;
		}

		// Prepare structures to send the old-new couples
		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> size;
		openfpm::vector<void *> ptr;

		// Array that will contain the couples, divided per processor
		openfpm::vector<openfpm::vector<size_t>> on_vs(vcl.getProcessingUnits());

		fillSendRecvStructs<size_t>(on_info, prc, size, ptr);

		// Send on_info
		vcl.sendrecvMultipleMessagesNBX(prc.size(), (size_t *)size.getPointer(), (size_t *)prc.getPointer(), (void **)ptr.getPointer(), on_receive, &on_vs, NONE);

		// Insert in the on_toup map the received couples
		for (size_t i = 0; i < vcl.getProcessingUnits(); i++)
		{
			if (i != vcl.getProcessUnitID() && on_vs.get(i).size() > 0) // redundant check 2nd arg in if
			{
				for (size_t j = 0; j < on_vs.get(i).size() - 1; j += 2) // -1 is useless
				{
					IdnProc nidpid = { on_vs.get(i).get(j + 1), i };
					on_toup.insert( { on_vs.get(i).get(j), nidpid });
				}
			}
		}

		// Update the glbi_map with the new ids and the processor info
		for (auto k : glbi_map)
		{
			auto search = on_toup.find(glbi_map.at(k.first).id); // fix with (k.second).id
			if (search != on_toup.end())
			{
				GlobalVInfo t = { (search->second).id, (search->second).pid };
				glbi_map.at(k.first) = t;
			}
		}

		// Vector of vertices global id I need info
		openfpm::vector<openfpm::vector<size_t>> vni(vcl.getProcessingUnits());

		// Map of re-mapping info
		std::unordered_map<size_t, size_t> rmi_m(vcl.getProcessingUnits()); // TODO elimin

		// Check which vertices I need to ask info about
		for (size_t i = 0; i < getNVertex(); ++i)
		{
			for (size_t j = 0; j < getNChilds(i); ++j)
			{
				// Here we get the global vertex id of all the children
				size_t vid = getChildDstGid(i, j);
				// We check which processor has the information about this vertex
				size_t pid = getInfoProc(vid);

				// if the vertex info is not in this processor I have to make a request to the right processor
				if (!vertexIsInThisGraph(vid) && pid != vcl.getProcessUnitID())
				{
					// add to requests
					vni.get(pid).add(vid);
				}
				else if (pid == vcl.getProcessUnitID())
				{
					// if info is in this processor add it in the map
					rmi_m.insert( { vid, glbi_map.at(vid).id });
				}
				else if (vertexIsInThisGraph(vid)) // check if it is needed, probably not because the glbi_map is update
				{
					// if the vertex is in this graph, add the new id in the map
					rmi_m.insert( { vid, glb2id.at(vid) });
				}
			}
		}

		// Array that will contain the requests from the other processors
		openfpm::vector<openfpm::vector<size_t>> req_rmi(vcl.getProcessingUnits());

		// Fill the structures for sendrecvMultipleMessagesNBX function
		fillSendRecvStructs<size_t>(vni, prc, size, ptr);

		// Send and receive requests
		vcl.sendrecvMultipleMessagesNBX(prc.size(), (size_t *)size.getPointer(), (size_t *)prc.getPointer(), (void **)ptr.getPointer(), on_receive, &req_rmi, NONE);

		// Re-mapping info map
		openfpm::vector<openfpm::vector<size_t>> rmi(vcl.getProcessingUnits());

		// Array that will contain the response to previous requests
		openfpm::vector<openfpm::vector<size_t>> resp_rmi(vcl.getProcessingUnits());

		// Prepare re-mapping info response
		for (size_t i = 0; i < req_rmi.size(); ++i)
		{
			for (size_t j = 0; j < req_rmi.get(i).size(); ++j)
			{
				resp_rmi.get(i).add(glbi_map.at(req_rmi.get(i).get(j)).id);
			}
		}

		// Fill the structures for sendrecvMultipleMessagesNBX function
		fillSendRecvStructs<size_t>(resp_rmi, prc, size, ptr);

		// Send responses
		vcl.sendrecvMultipleMessagesNBX(prc.size(), (size_t *)size.getPointer(), (size_t *)prc.getPointer(), (void **)ptr.getPointer(), on_receive, &rmi, NONE);

		// Add received info into re-mapping info map
		for (size_t i = 0; i < rmi.size(); ++i)
		{
			for (size_t j = 0; j < rmi.get(i).size(); ++j)
			{
				rmi_m.insert( { vni.get(i).get(j), rmi.get(i).get(j) });
			}
		}

		// Finally re-map the edges
		for (size_t i = 0; i < getNVertex(); ++i)
		{
			for (size_t s = 0; s < getNChilds(i); s++)
			{
				e_l.template get<e_map::vid>(i * v_slot + s) = rmi_m.at(getChildDstGid(i, s));
			}
		}
	}

	/*! \brief Initialize the fixed structure for global mapping. See glbiv for details.
	 *
	 */
	void initGlbimap()
	{
		size_t pid = vcl.getProcessUnitID();

		for (size_t i = 0; i < getNVertex(); ++i)
		{
			GlobalVInfo info;
			info.id = v_m.get(i).template get<v_info::id>();
			info.pid = pid;
			glbi_map.insert( { v_m.get(i).template get<v_info::gid>(), info });
		}
	}

	/*! \brief Get the processor id of the processor containing the vertex with global id vid
	 *
	 * \param vid vertex to know info about
	 * \return the processor id
	 */
	size_t getInfoProc(size_t vid)
	{
		for (size_t i = 0; i < fvtxdist.size() - 1; ++i)
		{
			if (vid >= (size_t)fvtxdist.get(i) && vid < (size_t)fvtxdist.get(i + 1))
			{
				return i;
			}
		}
		return vcl.getProcessingUnits() - 1;
	}

	/*! \brief Fill the prc, size and ptr structures with the data of vec
	 *
	 * \tparam type of the data inside vec
	 *
	 * \param vec vector of the data
	 * \param prc processor array of sendrecv function
	 * \param size size array of sendrecv function
	 * \param ptr pointers array of sendrecv function
	 */
	template<typename T>
	void fillSendRecvStructs(openfpm::vector<openfpm::vector<T>> &vec, openfpm::vector<size_t> &prc, openfpm::vector<size_t> &size, openfpm::vector<void *> &ptr)
	{
		// Reset sendrecv structures
		prc.clear();
		size.clear();
		ptr.clear();

		for (size_t i = 0; i < vec.size(); ++i)
		{
			if (vec.get(i).size() > 0 && i != vcl.getProcessUnitID())
			{
				prc.add(i);
				size.add(vec.get(i).size() * sizeof(T));
				ptr.add(vec.get(i).getPointer());
			}
		}
	}

public:

	//! Vertex typedef
	typedef V V_type;

	//! Edge typedef
	typedef E E_type;

	//! Object container for the vertex, for example can be encap<...> (map_grid or openfpm::vector)
	typedef typename openfpm::vector<V, Memory, layout_v, layout_v_base, grow_p, openfpm::vect_isel<V>::value>::container V_container;

	//! Object container for the edge, for example can be encap<...> (map_grid or openfpm::vector)
	typedef typename openfpm::vector<E, Memory, layout_e, layout_e_base, grow_p, openfpm::vect_isel<E>::value>::container E_container;

	/*! \brief It duplicate the graph
	 *
	 * \return a graph duplicate of the first
	 *
	 */

	DistGraph_CSR<V, E, Memory, layout_v, layout_e,layout_v_base,layout_e_base, grow_p> duplicate() const
	{
		DistGraph_CSR<V, E, Memory, layout_v, layout_e,layout_v_base,layout_e_base, grow_p> dup;

		dup.v_slot = v_slot;

		// duplicate all the structures

		dup.v.swap(v.duplicate());
		dup.v_m.swap(v_m.duplicate());
		dup.v_l.swap(v_l.duplicate());
		dup.glb2id = glb2id;
		dup.id2glb = id2glb;
		dup.glb2loc = glb2loc;
		dup.e.swap(e.duplicate());
		dup.e_m.swap(e_m.duplicate());
		dup.e_l.swap(e_l.duplicate());
		dup.e_invalid.swap(e_invalid.duplicate());
		dup.vtxdist.swap(vtxdist.duplicate());
		dup.fvtxdist.swap(fvtxdist.duplicate());
		return dup;
	}

	/*! \brief Constructor
	 *
	 * Constructor
	 *
	 */
	DistGraph_CSR(const DistGraph_CSR & dg) :
			vcl(create_vcluster())
	{
		this->operator=(dg);
	}

	/*! \brief Constructor
	 *
	 * Constructor
	 *
	 */
	DistGraph_CSR(DistGraph_CSR && dg)
	:vcl(create_vcluster())
	{
		this->operator=(dg);
	}

	/*! \brief Constructor
	 *
	 * Constructor
	 *
	 */
	DistGraph_CSR() :
			DistGraph_CSR(0, 16)
	{
	}

	/*! \brief Constructor
	 *
	 * Constructor
	 *
	 */
	DistGraph_CSR(size_t n_vertex) :
			DistGraph_CSR(n_vertex, 16)
	{
	}

	/*! \brief Constructor
	 *
	 * Constructor
	 *
	 */
	DistGraph_CSR(size_t n_vertex, size_t n_slot) :
			vcl(create_vcluster()), v_slot(n_slot)
	{
		// Creating n_vertex into the graph
		v.resize(n_vertex);
		// Creating n_vertex info objects into the graph
		v_m.resize(n_vertex);
		// Creating n_vertex adjacency list counters
		v_l.resize(n_vertex);
		// no edge set the counter to zero
		v_l.fill(0);
		// create one invalid edge
		e_invalid.resize(1);
		// init communication structures
		resetExchange();
	}

	/*! \brief Operator to access the decomposition vector
	 *
	 * \param v vector that will contain the decomposition
	 */
	void getDecompositionVector(openfpm::vector<idx_t> &v)
	{
		v.resize(vtxdist.size());

		for (size_t i = 0; i < vtxdist.size(); ++i)
		{
			v.get(i) = vtxdist.get(i);
		}
	}

	/*! \brief Operator to access the decomposition vector
	 *
	 * \return a pointer to the distribution vector
	 */
	openfpm::vector<idx_t>* getVtxdist()
	{
		return &vtxdist;
	}

	/*! \brief Operator to access the decomposition vector
	 *
	 * \param v vector that contains the decomposition
	 */
	void initDistributionVector(openfpm::vector<idx_t> &v)
	{
		vtxdist.resize(vcl.getProcessingUnits() + 1);
		fvtxdist.resize(vcl.getProcessingUnits() + 1);

		for (size_t i = 0; i < vtxdist.size(); ++i)
		{
			vtxdist.get(i) = v.get(i);
			fvtxdist.get(i) = v.get(i);
		}
	}

	/*! \brief Initialize the vtxdist and the fvtxdist
	 *
	 */
	void initDistributionVector()
	{
		updateVtxdist();

		fvtxdist.resize(vcl.getProcessingUnits() + 1);

		for (size_t i = 0; i < vtxdist.size(); ++i)
		{
			fvtxdist.get(i) = vtxdist.get(i);
		}
	}

	/*! \brief Copy constructor
	 *
	 * \param v_cl vcluster
	 * \param gg distributed graph to copy
	 *
	 */
	DistGraph_CSR(Vcluster & vcl, DistGraph_CSR<V, E, Memory> && g) :
			vcl(vcl)
	{
		swap(g);
	}

	/*! \brief Copy the graph
	 * 
	 * \param g distributed graph to copy
	 * 
	 * \return itself
	 *
	 */
	DistGraph_CSR<V, E, Memory> & operator=(DistGraph_CSR<V, E, Memory> && g)
	{
		swap(g);

		return *this;
	}

	/*! \brief Copy the graph
	 * 
	 * \param g graph to copy
	 * 
	 * \return itself
	 *
	 */
	DistGraph_CSR<V, E, Memory> & operator=(const DistGraph_CSR<V, E, Memory> & g)
	{
		swap(g.duplicate());

		return *this;
	}

	/*! \brief operator to access the vertex
	 *
	 * operator to access the vertex
	 *
	 * \tparam i property to access
	 * \param id of the vertex to access
	 *
	 * \return a reference to the vertex property
	 *
	 */
	template<unsigned int i> auto vertex_p(size_t id) -> decltype( v.template get<i>(id) )
	{
		return v.template get<i>(id);
	}

	/*! \brief Access the vertex
	 *
	 * \tparam i property to access
	 * \param id of the vertex to access
	 *
	 * \return a reference to the vertex property
	 *
	 */
	template<unsigned int i> auto vertex_p(grid_key_dx<1> id) -> decltype( v.template get<i>(id) )
	{
		return v.template get<i>(id);
	}

	/*! \brief Function to access the vertexes
	 *
	 * \param id of the vertex to access
	 *
	 * \return the vertex
	 *
	 */
	auto vertex(size_t id) -> decltype( v.get(id) )
	{
		return v.get(id);
	}

	/*! \brief operator to access the vertex
	 *
	 * operator to access the vertex
	 *
	 * \param id of the vertex to access
	 *
	 * \return the vertex
	 *
	 */
	auto vertex(grid_key_dx<1> id) -> decltype( v.get(id.get(0)) )
	{
		return v.get(id.get(0));
	}

	/*! \brief operator to access the vertex
	 *
	 * operator to access the vertex
	 *
	 * \param id of the vertex to access
	 *
	 * \return the vertex
	 *
	 */
	auto vertex(openfpm::vector_key_iterator id) -> decltype( v.get(0) )
	{
		return v.get(id.get());
	}

	/*! \brief Function to access the vertexes
	 *
	 * \param id of the vertex to access
	 *
	 * \return the vertex
	 *
	 */
	auto vertex(size_t id) const -> const decltype( v.get(id) )
	{
		return v.get(id);
	}

	/*! \brief operator to access the vertex
	 *
	 * operator to access the vertex
	 *
	 * \param id of the vertex to access
	 *
	 * \return the vertex
	 *
	 */
	auto vertex(grid_key_dx<1> id) const -> const decltype( v.get(id.get(0)) )
	{
		return v.get(id.get(0));
	}

	/*! \brief operator to access the vertex
	 *
	 * operator to access the vertex
	 *
	 * \param id of the vertex to access
	 *
	 * \return the vertex
	 *
	 */
	auto vertex(openfpm::vector_key_iterator id) const -> const decltype( v.get(0) )
	{
		return v.get(id.get());
	}

	/*! \brief operator to access the vertex info
	 *
	 * operator to access the vertex
	 *
	 * \param id of the vertex to access
	 *
	 * \return the vertex global id
	 *
	 */
	auto vertex_info(openfpm::vector_key_iterator id) const -> const decltype( v_m.get(0) )
	{
		return v_m.get(id.get());
	}

	/*! \brief Function to access the vertexes
	 *
	 * \param id GLOBAL id of the vertex to access
	 *
	 * \return the vertex
	 *
	 */
	auto getVertex(size_t id) -> decltype( v.get(id) )
	{

#ifdef SE_CLASS1

		if (glb2loc.find(id) == glb2loc.end())
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " The vertex with global id " << id << " is not in this sub-graph. Try to call reqVertex(" << id << ") and sync() first.\n";
			ACTION_ON_ERROR(DIST_GRAPH_ERROR);
		}

#endif

		return v.get(glb2loc.find(id)->second);
	}

	/*! \brief Function to access the vertexes
	 *
	 * \param id GLOBAL id of the vertex to access
	 *
	 * \return the vertex
	 *
	 */
	auto getVertex(size_t id) const -> const decltype( v.get(0) )
	{
		try
		{
			return v.get(glb2loc.at(id));
		} catch (const std::out_of_range& oor)
		{
			std::cerr << "The vertex with global id " << id << " is not in this sub-graph. Try to call reqVertex(" << id << ") and sync() first.\n";
		}

		return v.get(0);
	}

	/*! \brief operator to access the vertex position index by id property
	 *
	 * operator to access the vertex
	 *
	 * \param id id of the vertex to access
	 *
	 * \return id of the vertex
	 *
	 */
	size_t nodeById(size_t id) const
	{
		try
		{
			return glb2loc.at(id2glb.at(id));
		} catch (const std::out_of_range& oor)
		{
			std::cout << "Node not found by glb: " << id << std::endl;
		}

		return 0;
	}

	/*! /brief Get the first id of the graph
	 *
	 * \return the first id of the graph
	 */
	size_t firstId() const
	{
		return vtxdist.get(vcl.getProcessUnitID());
	}

	/*! /brief Get the last id of the graph
	 *
	 * \return the last id of the graph
	 */
	size_t lastId() const
	{
		return vtxdist.get(vcl.getProcessUnitID() + 1) - 1;
	}

	/*! \brief Get the id of a vertex given its index position
	 *
	 * \param i position of the vertex
	 * \return the id of the vertex
	 */
	size_t getVertexId(size_t i) const
	{
		return v_m.get(i).template get<v_info::id>();
	}

	/*! \brief Get the id of a vertex given its index position
	 *
	 * \param i position of the vertex
	 * \return the id of the vertex
	 */
	size_t getVertexGlobalId(size_t i) const
	{
		return v_m.get(i).template get<v_info::gid>();
	}

	/*! \brief Check if the vertex with GLOBAL id is in this graph
	 *
	 * \param id global id of the vertex
	 * \return true if vertex with gid is in this graph, false otherwise
	 */
	bool vertexIsInThisGraph(size_t id)
	{
		auto search = glb2id.find(id);

		if (search != glb2id.end())
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/*! \brief operator to update all the hashmap
	 *
	 * \param n new vertex id
	 * \param g global vertex id
	 * \param l local vertex id
	 *
	 * \tparam i id of the property storing the id
	 *
	 */
	void map_v(size_t n, size_t g, size_t l)
	{
		id2glb.insert( { n, g });
		glb2id.insert( { g, n });
		glb2loc.insert( { g, l });
		v_m.get(l).template get<v_info::id>() = n;
	}

	/*! \brief operator to clear the whole graph
	 *
	 * operator to clear all
	 *
	 */
	void clear()
	{
		v.clear();
		v_m.clear();
		e.clear();
		id2glb.clear();
		glb2id.clear();
		glb2loc.clear();
		v_l.clear();
		e_l.clear();
		e_invalid.clear();
	}

	/*! \brief Access the edge
	 *
	 * \tparam i property to access
	 * \param id of the edge to access
	 *
	 * \return a reference to the edge property
	 *
	 */
	template<unsigned int i> auto edge_p(grid_key_dx<1> id) -> decltype ( e.template get<i>(id) )
	{
		return e.template get<i>(id);
	}

	/*! \brief Access the edge
	 *
	 * \tparam i property to access
	 * \param id of the edge to access
	 *
	 * \return a reference to the edge property
	 *
	 */
	template<unsigned int i> auto edge_p(size_t id) -> decltype ( e.template get<i>(id) )
	{
		return e.template get<i>(id);
	}

	/*! \brief Access the edge
	 *
	 * \param id of the edge to access
	 *
	 * \return a reference to the edge
	 *
	 */
	auto edge(grid_key_dx<1> id) const -> const decltype ( e.get(id.get(0)) )
	{
		return e.get(id.get(0));
	}

	/*! \brief operator to access the edge
	 *
	 * \param ek key of the edge
	 *
	 * \return a reference to the edge
	 *
	 */
	auto edge(edge_key ek) const -> const decltype ( e.get(0) )
	{
		return e.get(e_l.template get<e_map::eid>(ek.pos * v_slot + ek.pos_e));
	}

	/*! \brief operator to access the edge
	 *
	 * \param ek key of the edge
	 *
	 * \param return a reference to the edge
	 *
	 */
	auto getEdge(edge_key ek) const -> const decltype ( e.get(0) )
	{
		size_t v;

		try
		{
			v = glb2loc.at(ek.pos);
		} catch (const std::out_of_range& oor)
		{
			std::cout << "The source vertex of this edge is not in this graph.\n";
		}

		return e.get(e_l.template get<e_map::eid>(v * v_slot + ek.pos_e));
	}

	/*! \brief operator to access the edge
	 *
	 * operator to access the edge
	 *
	 * \param id of the edge to access
	 *
	 * \return a reference to the edge
	 *
	 */
	auto edge(size_t id) const -> const decltype ( e.get(id) )
	{
		return e.get(id);
	}

	/*! \brief Return the number of children of a vertex
	 *
	 * \param c Child id
	 *
	 * \return the number of children
	 *
	 */
	inline size_t getNChilds(size_t c) const
	{
		return v_l.template get<0>(c);
	}

	/*! \brief Return the number of childs of a vertex
	 *
	 * \param c child id
	 *
	 * \return the number of childs
	 *
	 */
	inline size_t getNChilds(typename openfpm::vector<V, Memory, layout_v, layout_v_base, grow_p, openfpm::vect_isel<V>::value>::iterator_key & c)
	{
		return v_l.template get<0>(c.get());
	}

	/*! \brief Return the number of children of a vertex given its global id
	 *
	 * \param v vertex global id
	 *
	 * \return the number of children
	 *
	 */
	inline size_t getNEdge(size_t v) const
	{
		try
		{
			v = glb2loc.at(v);
		} catch (const std::out_of_range& oor)
		{
			std::cout << "The source vertex of this edge is not in this graph.\n";
		}

		return v_l.template get<0>(v);
	}

	/*! \brief Get the vertex edge
	 *
	 * \param v vertex
	 * \param v_e edge id
	 *
	 * \return the edge
	 *
	 */
	inline auto getChildEdge(size_t v, size_t v_e) -> decltype(e.get(0))
	{
		return e.get(e_l.template get<e_map::eid>(v * v_slot + v_e));
	}

	/*! \brief Get the vertex edge info
	 *
	 * \param v vertex
	 * \param v_e edge id
	 *
	 * \return the id of the edge
	 *
	 */
	inline auto getChildInfo(size_t v, size_t v_e) -> decltype(e_m.get(0))
	{
		return e_m.get(e_l.template get<e_map::eid>(v * v_slot + v_e));
	}

	/*! \brief Get the vertex edge given the vertex global id as source
	 *
	 * \param v vertex global id
	 * \param v_e edge id
	 *
	 * \return the edge
	 *
	 */
	inline auto getEdge(size_t v, size_t v_e) -> decltype(e.get(0))
	{
		try
		{
			v = glb2loc.at(v);
		} catch (const std::out_of_range& oor)
		{
			std::cout << "The source vertex of this edge is not in this graph.\n";
		}

		return e.get(e_l.template get<e_map::eid>(v * v_slot + v_e));
	}

	/*! \brief Get the child edge
	 *
	 * \param v node
	 * \param i child at position i
	 *
	 * \return the edge id that connect v with the target at position i
	 *
	 */
	inline size_t getChild(size_t v, size_t i) const
	{
#ifdef SE_CLASS1
		if (i >= v_l.template get<0>(v))
		{
			std::cerr << "Error " << __FILE__ << " line: " << __LINE__ << "    vertex " << v << " does not have edge " << i << " on processor " << vcl.getProcessUnitID() << "\n";
		}

		if (e.size() <= e_l.template get<e_map::eid>(v * v_slot + i))
		{
			std::cerr << "Error " << __FILE__ << " " << __LINE__ << " vertex " << v << " does not have edge " << i << " on processor " << vcl.getProcessUnitID() << " (" << e.size() << "<=" << e_l.template get<e_map::eid>(v * v_slot + i) << ")\n";
		}
#endif
		// Get the edge id
		return e_l.template get<e_map::vid>(v * v_slot + i);
	}

	/*! \brief Get the child edge
	 *
	 * \param i id of the child
	 *
	 * \return the edge id that connect v with the target at position i
	 *
	 */
	inline size_t getChild(size_t i) const
	{
		// Get the edge id
		return e_l.template get<e_map::vid>(i);
	}

	/*! \brief Get the child edge
	 *
	 * \param v node
	 * \param i child at position i
	 *
	 * \return the target i connected by an edge node, for the node v
	 *
	 */
	inline size_t getChild(typename openfpm::vector<V, Memory, layout_v, layout_v_base, grow_p, openfpm::vect_isel<V>::value>::iterator_key & v, size_t i)
	{
#ifdef DEBUG
		if (i >= v_l.template get<0>(v.get()))
		{
			std::cerr << "Error " << __FILE__ << " line: " << __LINE__ << "    vertex " << v.get() << " does not have edge " << i << "\n";
		}

		if (e.size() <= e_l.template get<e_map::eid>(v.get() * v_slot + i))
		{
			std::cerr << "Error " << __FILE__ << " " << __LINE__ << " vertex " << v.get() << " does not have edge " << i << "\n";
		}
#endif

		// Get the edge id
		return e_l.template get<e_map::vid>(v.get() * v_slot + i);
	}

	/*! \brief Add vertex vrt with global id and id properties
	 *
	 * \param vrt vertex object to add
	 * \param gid global id, unique in global graph
	 * \param id id, unique n global graph
	 */
	inline void add_vertex(const V & vrt, size_t id, size_t gid)
	{
		// Create vertex info object
		v_info vm;
		vm.template get<v_info::id>() = id;
		vm.template get<v_info::gid>() = gid;

		// Add the vertex info
		v_m.add(vm);

		// Add the vertex
		v.add(vrt);

		// Update id to global map
		id2glb.insert( { id, gid });

		// Update global id to local index
		glb2loc.insert( { gid, v.size() - 1 });

		// Update global id to id
		glb2id.insert( { gid, id });

		// Set the number of adjacent vertex for this vertex to 0
		v_l.add(0ul);

		// Add a slot for the vertex adjacency list
		e_l.resize(e_l.size() + v_slot);
	}

	/*! \brief Add vertex vrt with global id and id properties
	 *
	 * \param vrt vertex object to add
	 * \param id of the vertex
	 * \param gid global id, unique in global graph
	 */
	template<unsigned int dim, typename Mem> inline void add_vertex(const encapc<dim, V, Mem> & vrt, size_t id, size_t gid)
	{
		// Create vertex info object
		v_info vm;
		vm.template get<v_info::id>() = id;
		vm.template get<v_info::gid>() = gid;

		// Add the vertex info
		v_m.add(vm);

		// Add the vertex
		v.add(vrt);

		// Update id to global map
		id2glb.insert( { id, gid });

		// Update global id to local index
		glb2loc.insert( { gid, v.size() - 1 });

		// Update global id to id
		glb2id.insert( { gid, id });

		// Set the number of adjacent vertex for this vertex to 0
		v_l.add(0ul);

		// Add a slot for the vertex adjacency list
		e_l.resize(e_l.size() + v_slot);
	}

	/*! \brief Add vertex vrt with global id and id properties
	 *
	 * \param vrt vertex object to add
	 * \param gid global id, unique in global graph
	 */
	inline void add_vertex(const V & vrt, size_t gid)
	{
		add_vertex(vrt, gid, gid);
	}

	/*! \brief map global id to local id
	 *
	 * \param g global id
	 * \param l local index
	 * \param i id
	 */
	void setGlobalMap(size_t g, size_t l, size_t i)
	{
		v_m.template get<v_info::id>(l) = i;

		v_m.template get<v_info::gid>(l) = g;

		// Set global id to local index
		glb2loc.insert( { g, l });

		// Set global id to id
		glb2id.insert( { g, i });

		// Set id to global map
		id2glb.insert( { i, g });
	}

	inline auto addEdge(size_t v1, size_t v2, size_t srcgid, size_t dstgid) -> decltype(e.get(0))
	{
		// add an edge
		long int id_x_end = addEdge_<NoCheck>(v1, v2);
		// If there is not edge return an invalid edge, is a kind of stub object
		if (id_x_end == NO_EDGE)
			return e_invalid.get(0);

		// set source and destination global ids of the edge
		e_m.template get<e_info::sgid>(id_x_end) = srcgid;
		e_m.template get<e_info::dgid>(id_x_end) = dstgid;

		// return the edge to change the properties
		return e.get(id_x_end);
	}

	inline auto addEdge(size_t v1, size_t v2, size_t srcgid, size_t dstgid, const E & ed) -> decltype(e.get(0))
	{
		// add an edge
		long int id_x_end = addEdge_<NoCheck>(v1, v2);
		// If there is not edge return an invalid edge, is a kind of stub object
		if (id_x_end == NO_EDGE)
			return e_invalid.get(0);

		// add in e_l the edge properties
		e.set(id_x_end, ed);

		// set source and destination global ids of the edge
		e_m.template get<e_info::sgid>(id_x_end) = srcgid;
		e_m.template get<e_info::dgid>(id_x_end) = dstgid;

		// return the edge to change the properties
		return e.get(id_x_end);
	}

	template<unsigned int dim, typename Mem, typename Mem1> inline auto addEdge(size_t v1, size_t v2, const encapc<dim, E, Mem> & ed, const encapc<dim, e_info, Mem1> & ei) -> decltype(e.get(0))
	{
		// add an edge
		long int id_x_end = addEdge_<NoCheck>(v1, v2);
		// If there is not edge return an invalid edge, is a kind of stub object
		if (id_x_end == NO_EDGE)
			return e_invalid.get(0);

		// set the edge object and the edge info object
		e.set(id_x_end, ed);
		e_m.set(id_x_end, ei);

		// return the edge to change the properties
		return e.get(id_x_end);
	}

	inline auto addEdge(size_t v1, size_t v2, const E & ed, const e_info & ei) -> decltype(e.get(0))
	{
		// add an edge
		long int id_x_end = addEdge_<NoCheck>(v1, v2);
		// If there is not edge return an invalid edge, is a kind of stub object
		if (id_x_end == NO_EDGE)
			return e_invalid.get(0);

		// set the edge object and the edge info object
		e.set(id_x_end, ed);
		e_m.set(id_x_end, ei);

		// return the edge to change the properties
		return e.get(id_x_end);
	}

	/*! \brief Get the global id of edge's source vertex
	 *
	 * \param v1 source vertex of the edge
	 * \param s n child of vertex v1
	 * \return global id of the source vertex
	 */
	size_t getChildSrcGid(size_t v1, size_t s)
	{
		size_t eid = e_l.template get<e_map::eid>(v1 * v_slot + s);
		return e_m.template get<e_info::sgid>(eid);
	}

	/*! \brief Get the global id of edge's destination vertex
	 *
	 * \param v1 source vertex of the edge
	 * \param s n child of vertex v1
	 * \return global id of the destination vertex
	 */
	size_t getChildDstGid(size_t v1, size_t s)
	{
		size_t eid = e_l.template get<e_map::eid>(v1 * v_slot + s);
		return e_m.template get<e_info::dgid>(eid);
	}

	/*! \brief Add an edge between vertices v1 end v2, needs syncEdge() to complete the action
	 *
	 * \param v1 source vertex of the edge
	 * \param v2 destination vertex of the edge
	 */
	template<typename CheckPolicy = NoCheck> inline void add_edge(size_t v1, size_t v2)
	{
		//if the source vertex is not in this graph, this processor doesn't need to do anything
		if (!vertexIsInThisGraph(v1))
		{
			return;
		}

		//if the destination vertex is not in this graph, this processor has to request it and add the edge to a queue
		if (!vertexIsInThisGraph(v2))
		{
			reqVertex(v2);
			EdgeReq er = { v1, v2, 0, 0 };
			e_queue.add(er);

			return;
		}

		// add an edge
		long int id_x_end = addEdge_<CheckPolicy>(glb2loc.at(v1), glb2id.at(v2));

		// If there is not edge return an invalid edge, is a kind of stub object
		if (id_x_end == NO_EDGE)
			return;

		// set source and destination ids of the edge
		e_m.get(id_x_end).template get<e_info::sgid>() = v1;
		e_m.get(id_x_end).template get<e_info::dgid>() = v2;
	}

	/*! \brief Execute a synchronization through processor to finalize the add of the edges requested in the e_queue
	 *
	 */
	void syncEdge()
	{
		// retrieve ghosts necessary for edge adding
		sync();

		// update with real values of edges, we don't add the edge here because it will modify edge's array length
		for (size_t i = 0; i < e_queue.size(); ++i)
		{
			e_queue.get(i).v1n = glb2loc.at(e_queue.get(i).v1);
			e_queue.get(i).v2n = glb2id.at(e_queue.get(i).v2);
		}

		deleteGhosts();

		for (auto req : e_queue)
		{
			// add an edge
			long int id_x_end = addEdge_<>(req.v1n, req.v2n);

			// If there is not edge return an invalid edge, is a kind of stub object
			if (id_x_end == NO_EDGE)
				return;

			// set source and destination ids of the edge
			e_m.get(id_x_end).template get<e_info::sgid>() = req.v1;
			e_m.get(id_x_end).template get<e_info::dgid>() = req.v2;
		}

		e_queue.clear();

	}

	/*! \brief Swap the memory of g with this graph
	 *
	 * Swap the memory of g with this graph, it is basically used
	 * for move semantic
	 *
	 * \param g The source graph
	 */
	inline void swap(DistGraph_CSR<V, E> & g)
	{
		// switch the memory
		v.swap(g.v);
		v_m.swap(g.v_m);
		e.swap(g.e);
		e_m.swap(g.e_m);
		v_l.swap(g.v_l);
		glb2id.swap(g.glb2id);
		id2glb.swap(g.id2glb);
		glb2loc.swap(g.glb2loc);
		e_l.swap(g.e_l);
		e_invalid.swap(g.e_invalid);
		vtxdist.swap(g.vtxdist);
		fvtxdist.swap(g.fvtxdist);

		size_t v_slot_tmp = v_slot;
		v_slot = g.v_slot;
		g.v_slot = v_slot_tmp;
	}

	/*! \brief Swap the memory of g with this graph
	 *
	 * Swap the memory of g with this graph, it is basically used
	 * for move semantic
	 *
	 * \param g The source graph
	 */
	inline void swap(DistGraph_CSR<V, E> && g)
	{
		// switch the memory
		v.swap(g.v);
		v_m.swap(g.v_m);
		e.swap(g.e);
		e_m.swap(g.e_m);
		v_l.swap(g.v_l);
		glb2id.swap(g.glb2id);
		id2glb.swap(g.id2glb);
		glb2loc.swap(g.glb2loc);
		e_l.swap(g.e_l);
		e_invalid.swap(g.e_invalid);
		vtxdist.swap(g.vtxdist);
		fvtxdist.swap(g.fvtxdist);

		size_t v_slot_tmp = v_slot;
		v_slot = g.v_slot;
		g.v_slot = v_slot_tmp;
	}

	/*! \brief Get the vertex iterator
	 *
	 * Get the vertex iterator
	 *
	 * \return an iterator to iterate through all the vertex
	 *
	 */
	inline auto getVertexIterator() const -> decltype(v.getIterator())
	{
		return v.getIterator();
	}

	/*! \brief Get the vertex iterator
	 *
	 * Get the vertex iterator
	 *
	 * \return an iterator to iterate through all the edges
	 *
	 */
	inline edge_iterator<DistGraph_CSR<V, E, Memory>> getEdgeIterator() const
	{
		return edge_iterator<DistGraph_CSR<V, E, Memory>>(*this);
	}

	/*! \brief Return the number of the vertices in this subgraph
	 *
	 * \return the number of vertex
	 *
	 */
	inline size_t getNVertex() const
	{
		return v.size();
	}

	/*! \brief Return the total number of the vertices
	 *
	 * \return the number of vertex
	 *
	 */
	inline size_t getTotNVertex() const
	{
		return vtxdist.get(vcl.getProcessingUnits());
	}

	/*! \brief Return the number of edges
	 *
	 * \return the number of edges
	 *
	 */
	inline size_t getNEdge() const
	{
		return e.size();
	}

	/*! \brief Once added all the vertices this function must be called to initialize all the properties, useless if a graph factory is used
	 *
	 */
	void init()
	{
		initDistributionVector();

		initGlbimap();

		remap();
	}

	/*! \brief Check if a vertex is a ghost vertex (not belonging to this processor)
	 *
	 * \param id id of the vertex to check
	 * \return true if it is a ghost vertex
	 */
	bool isGhost(size_t id)
	{
		auto search = ghs_map.find(id);

		if (search != ghs_map.end())
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/*! \brief Remove all the ghosts from this graph
	 *
	 */
	void deleteGhosts()
	{
		if (ghs_map.size() == 0)
			return;

		// Prepare new graph
		DistGraph_CSR recv_g;

		size_t fpos = getNVertex() - ghs_map.size();
		size_t epos = 0;

		// Reset global to local map
		for (auto gh : ghs_map)
		{
			epos += getNChilds(glb2loc.at(gh.first));
			id2glb.erase(glb2id.at(gh.first));
			glb2loc.erase(gh.first);
			glb2id.erase(gh.first);

		}

		//resize all structures to delete the ghosts
		v.resize(fpos);
		v_m.resize(fpos);
		v_l.resize(fpos);
		e.resize(e.size() - epos);
		e_m.resize(e_m.size() - epos);
		e_l.resize(fpos * v_slot);

		// Clear ghosts map
		ghs_map.clear();
	}

	/*! \brief Prepare to send vertex i from the local processor to the target processor
	 *
	 * \tparam toRemove if true the vertex will be deleted from this graph once it has been sent
	 *
	 * \param i vertex to send
	 * \param t target processor
	 */
	template<bool toRemove = true> //TODO make it private and create wrapper in public
	void q_move(size_t i, size_t t)
	{
		// Check if a 'useless' move has been requested
		if (t == vcl.getProcessUnitID())
		{
			std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " target processor is equal the local processor\n";
			return;
		}

		// If the array for t processor is empty, set it to not empty
		if (sgp.get(t).isEmpty)
			sgp.get(t).isEmpty = false;

		// Add the vertex to the vertex send buffer
		sgp.get(t).send_v.add(vertex(i));

		// Add the vertex to the vertex send buffer
		sgp.get(t).send_v_m.add(v_m.get(i));

		// Add the info about the number of children
		sgp.get(t).send_es.add(getNChilds(i));

		for (size_t s = 0; s < getNChilds(i); s++)
		{
			// Add edge value and object to the respective arrays
			sgp.get(t).send_e.add(getChildEdge(i, s));
			sgp.get(t).send_e_m.add(getChildInfo(i, s));
			sgp.get(t).send_el.add(getChild(i, s));
		}

		// If the vertex has to be removed after the send add its index id to the v_td array
		if (toRemove)
			v_td.add(i);
	}

	/*! \brief Check if the move queue is empty
	 *
	 * \return true if the move queue is empty
	 */
	bool moveQueueIsEmpty()
	{
		size_t exs = 0;
		for (size_t i = 0; i < vcl.getProcessingUnits(); ++i)
			if (sgp.get(i).isEmpty)
				exs++;

		if (exs == 0)
			return true;

		return false;
	}

	/*! \brief Redistribute function that wraps different stages of the redistribution
	 *
	 */
	void redistribute()
	{
		if (glbi_map.size() == 0)
			initGlbimap();

		deleteGhosts();

		exchangeVertices<false>();

		updateVtxdist();

		remap();
	}

	/*! \brief Put a vertex request in queue
	 *
	 * \param gid global id of the vertex to request
	 */
	void reqVertex(size_t gid)
	{
		// if not initialized, prepare vr_queue
		if (vr_queue.size() == 0)
			vr_queue.resize(vcl.getProcessingUnits());

		if (!vertexIsInThisGraph(gid))
		{
			vr_queue.get(getInfoProc(gid)).add(gid);
		}

	}

	/*! \brief Execute all vertex requests and add them as ghosts inside this graph, they will be available until a redistribution is executed
	 *
	 */
	void sync()
	{
		// If not initialized, prepare global informations map
		if (glbi_map.size() == 0)
			initGlbimap();

		// if not initialized, prepare vr_queue
		if (vr_queue.size() == 0) // TODO remove check on size
			vr_queue.resize(vcl.getProcessingUnits());

		// Arrays that will contain temporary requests and responses during communications TODO move to global private and reset method
		openfpm::vector<openfpm::vector<size_t>> resp(vcl.getProcessingUnits());
		openfpm::vector<openfpm::vector<size_t>> reqs(vcl.getProcessingUnits());

		// Prepare structures for communication
		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> size;
		openfpm::vector<void *> ptr;

		// Fill the structures for sendrecvMultipleMessagesNBX function
		fillSendRecvStructs<size_t>(vr_queue, prc, size, ptr);

		// Send/receive requests for info about needed vertices
		vcl.sendrecvMultipleMessagesNBX(prc.size(), (size_t *)size.getPointer(), (size_t *)prc.getPointer(), (void **)ptr.getPointer(), on_receive, &resp, NONE);

		// Prepare responses with the containing processors of requested vertices
		for (size_t i = 0; i < resp.size(); ++i)
		{
			reqs.get(i).clear();

			for (size_t j = 0; j < resp.get(i).size(); ++j)
			{
				try
				{
					reqs.get(i).add(glbi_map.at(resp.get(i).get(j)).pid);
				} catch (const std::out_of_range& oor)
				{
					std::cout << resp.get(i).get(j) << " not found in global info map (proc: " << vcl.getProcessUnitID() << ")\n";
				}
			}
		}

		// Fill the structures for sendrecvMultipleMessagesNBX function
		fillSendRecvStructs<size_t>(reqs, prc, size, ptr);

		// Reset response array TODO clear and resize not needed
		resp.clear();
		resp.resize(vcl.getProcessingUnits());

		// Send/receive responses with the containing processors of requested vertices
		vcl.sendrecvMultipleMessagesNBX(prc.size(), (size_t *)size.getPointer(), (size_t *)prc.getPointer(), (void **)ptr.getPointer(), on_receive, &resp, NONE);

		// Clear requests array
		reqs.clear();
		reqs.resize(vcl.getProcessingUnits());

		// Prepare vertices requests
		for (size_t i = 0; i < vcl.getProcessingUnits(); ++i)
		{
			// If the info has been retrieved from other processors take it from resp
			if (i != vcl.getProcessUnitID())
			{
				for (size_t j = 0; j < vr_queue.get(i).size(); ++j)
				{
					reqs.get(resp.get(i).get(j)).add(vr_queue.get(i).get(j));
				}
			}
			// Otherwise take it from the global info map of this processor
			else
			{
				for (size_t j = 0; j < vr_queue.get(i).size(); ++j)
				{
					reqs.get(glbi_map.at(vr_queue.get(i).get(j)).pid).add(vr_queue.get(i).get(j));
				}
			}
		}

		// Fill the structures for sendrecvMultipleMessagesNBX function
		fillSendRecvStructs<size_t>(reqs, prc, size, ptr);

		// Reset response array
		resp.clear();
		resp.resize(vcl.getProcessingUnits());

		// Send/receive vertices requests
		vcl.sendrecvMultipleMessagesNBX(prc.size(), (size_t *)size.getPointer(), (size_t *)prc.getPointer(), (void **)ptr.getPointer(), on_receive, &resp, NONE);

		for (size_t i = 0; i < resp.size(); ++i)
		{
			for (size_t j = 0; j < resp.get(i).size(); ++j)
			{
				try
				{
					q_move<false>(glb2loc.at(resp.get(i).get(j)), i);
				} catch (const std::out_of_range& oor)
				{
					std::cout << resp.get(i).get(j) << " not found in global to local map (proc: " << vcl.getProcessUnitID() << ")\n";
				}
			}
		}

		// Send and receive vertices, the received ones will be added to the graph as ghosts
		exchangeVertices<true>();

		// Empty the requests list
		vr_queue.clear();
	}
};

/*! \brief Simplified implementation of DistGraph_CSR
 *
 * Used when DistGraph_CSR is used as a default template argument to avoid 7 arguments
 *
 * [Example]
 *
 * template<template<typename,typename> class T=DistGraph_CSR_s>
 * class cool_structure
 * {
 * 		T<Vertex,Edge> graph
 * }
 *
 * only 2 parameter are needed, if you use DistGraph_CSR you have to define 7 regardless that
 * DistGraph_CSR has some default template
 *
 * template<template<typename,typename> class T=DistGraph_CSR>
 * class cool_structure
 * {
 * 		T<Vertex,Edge> graph
 * }
 *
 * THIS DO NOT COMPILE
 *
 */
template<typename V, typename E>
class DistGraph_CSR_s: public DistGraph_CSR<V, E>
{

};

#endif /* DIST_MAP_GRAPH_HPP_ */
