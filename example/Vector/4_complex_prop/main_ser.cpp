/*!
 *
 * \page Vector_4_complex_prop_ser Vector 4 property serialization
 *
 *
 * [TOC]
 *
 *
 * # Vector 4 property serialization # {#vector_example_cp_ser}
 *
 *
 * This example show how we can use complex properties in a vector when we use a custom structure with
 * a pointer inside. In particular we will show how to construct the serialization and de-serialization
 * methods for the structure my_struct defined here
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp my structure
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp my structure end
 *
 */

#include "Vector/vector_dist.hpp"

//! \cond [my structure] \endcond

struct my_struct
{
	//! C string size
	size_t size;

	//! C string pointer
	char * ptr;

	//! C++ string
	std::string str;

	//! vector
	openfpm::vector<int> v;

//! \cond [my structure] \endcond

public:

	//! Functions to check if the packing object is complex
	static bool pack()	{return true;}


	//! \cond [con and dest] \endcond

	my_struct()
	{
		size = 0;
		ptr = NULL;
	}

	~my_struct()
	{
		if (ptr != NULL)
			delete [] ptr;
	}

	//! \cond [con and dest] \endcond

	//! \cond [pack request] \endcond

	//! Serialization request
	template<int ... prp> inline void packRequest(size_t & req) const
	{
		req += 2*sizeof(size_t) + size + str.size()+1;
		v.packRequest(req);
	}

	//! \cond [pack request] \endcond

	//! \cond [pack serialize] \endcond

	//! Serialize the data structure
	template<typename Memory, int ... prp> inline void pack(ExtPreAlloc<Memory> & mem, Pack_stat & sts) const
	{
		//! \cond [packer ser pack] \endcond

		//Serialize the number that determine the C-string size
		Packer<size_t, Memory>::pack(mem,size,sts);

		//! \cond [packer ser pack] \endcond

		//! \cond [pack ser copy] \endcond

		// Allocate and copy the string
		mem.allocate(size);
		void * t_ptr = mem.getPointer();
		memcpy(t_ptr,ptr,size);

		//! \cond [pack ser copy] \endcond

		//! \cond [pack ser other] \endcond

		//Serialize the number that determine the C++-string str.size()+1 given by the null terminator
		Packer<size_t, Memory>::pack(mem,str.size()+1,sts);

		// Allocate and copy the string
		mem.allocate(str.size()+1);
		char * t_ptr2 = (char *)mem.getPointer();
		memcpy(t_ptr2,str.c_str(),str.size());

		// Add null terminator
		t_ptr2[str.size()] = 0;

		Packer<decltype(v), Memory>::pack(mem,v,sts);

		//! \cond [pack ser other] \endcond
	}

	//! \cond [pack serialize] \endcond

	//! \cond [unpack de-serialize] \endcond

	//! De-serialize the data structure
	template<typename Memory, int ... prp> inline void unpack(ExtPreAlloc<Memory> & mem, Unpack_stat & ps)
	{
		//! \cond [unpacker ser pack] \endcond

		// unpack the size of the C string
		Unpacker<size_t, Memory>::unpack(mem,size,ps);

		//! \cond [unpacker ser pack] \endcond

		//! \cond [unpacker ser create] \endcond

		// Allocate the C string
		ptr = new char[size];

		// source pointer
		char * ptr_src = (char *)mem.getPointerOffset(ps.getOffset());

		// copy from the buffer to the destination
		memcpy(ptr,ptr_src,size);

		// update the pointer
		ps.addOffset(size);

		//! \cond [unpacker ser create] \endcond

		//! \cond [unpacker ser other] \endcond

		// get the the C++ string size
		size_t cpp_size;
		Unpacker<size_t, Memory>::unpack(mem,cpp_size,ps);

		// Create the string from the serialized data (de-serialize)
		char * ptr_src2 = (char *)mem.getPointerOffset(ps.getOffset());
		str = std::string(ptr_src2);
		ps.addOffset(cpp_size);

		// Unpack the vector
		Unpacker<decltype(v), Memory>::unpack(mem,v,ps);

		//! \cond [unpacker ser other] \endcond
	}

	//! \cond [unpack de-serialize] \endcond

//! \cond [my structure end] \endcond

};

//! \cond [my structure end] \endcond

/*!
 *
 * \page Vector_4_complex_prop_ser Vector 4 property serialization
 *
 * In order to make my_struct serializable we need 3 methods
 *
 * * **packRequest** This method indicate how many byte are needed to serialize this structure
 * * **pack** This method serialize the data-structure
 * * **unpack** This method de-serialize the data structure
 *
 * ### packRequest ###
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp pack request
 *
 * This function calculate the size in byte to serialize this structure in this case we serialize
 * 2 numbers so 2*sizeof(size_t). One C string of size "size" and one C++ string of size str.size(). Because
 * std::string has a constructor from null terminating string we add the null terminator to easy construct an
 * std::string. This mean that
 * the C++ string will be serialized in str.size()+1 bytes. The result must be summed to counter req that is used
 * to allocate a buffer big enough for serialization. The function is template and has a variadic argument
 * "int ... prp" this can be ignored
 *
 * ### pack ###
 *
 * Here we serialize the object. The function give as argument
 *
 * * **ExtPreAlloc<Memory>** mem The memory where we have to serialize the information
 * * **Pack_stat** unused
 * * **int ... prp** unused
 *
 * The only important parameter is the object mem where we will serialize the information
 *
 * We first pack the number that indicate the size of the C string. A convenient way
 * is use the function
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp packer ser pack
 *
 * The function **Packer<decltype(object),Memory>::pack(object,sts)** work if **object** is a
 * fundamental C++ type, a struct that does not contain pointers, any object that
 * implement the interface pack,unpack,packRequest (like my_struct). Most openfpm
 * data-structure like openfpm::vector and grid implement such interface and
 * can be used directly with **Packer<...>::pack(...)**.
 * After packed the size of the string we have to serialize the content of the pointer string.
 * unfortunately there is no Packer<...>::pack(...) function to do this and must be
 * manually done. This can be done with the following steps
 *
 * * Allocate the memory to copy the value of the string
 * * Get the pointer to the memory allocated
 * * copy the content of the string into the allocated memory
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp pack ser copy
 *
 * For the C++ string the concept is the same, the difference is that we put a null terminator at the end
 * of the serialized string
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp pack ser other
 *
 * ### unpack ###
 *
 * Here we de-serialize the object. The function take a reference to object
 *
 * * **ExtPreAlloc<Memory>** mem, contain the serialized information that must be de-serialized
 * * **Pack_stat** contain the offset to the memory to de-serialize
 * * **int ... prp** unused
 *
 * De-serialization must be in general done in the same order as we serialized
 * We first unpack the number that indicate the size of the C string. A convenient way
 * is to use the function
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp unpacker ser pack
 *
 * the variable **size** contain now the size of the packed string we can now
 *
 * * create memory with **new** that store the C string
 * * Get the pointer to the serialized information
 * * copy the string from mem into the created memory
 * * update the offset pointer
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp unpacker ser create
 *
 * For the C++ string is just the same
 *
 * * We get the size of the string
 * * we create an std::string out of the null terminating serialized string
 * * we assign the created std::string to str
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp unpacker ser other
 *
 * ### Constructor and destructor ###
 *
 * Constructor and destructor are not releated to serialization and de-serialization concept.
 * But on how my_struct is constructed and destructed in order to avoid memory
 * leak/corruption. In the constructor we set ptr to NULL in the destructor we
 *  destroy the pointer (if different from NULL) to avoid memory leak
 *
 * \snippet Vector/4_complex_prop/main_ser.cpp con and dest
 *
 */

int main(int argc, char* argv[])
{
	/*!
	 *
	 * \page Vector_4_complex_prop_ser Vector 4 property serialization
	 *
	 *
	 * ## Initialization and vector creation ##
	 *
	 * After we initialize the library we can create a vector with complex properties
	 * with the following line
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vect create
	 *
	 * In this this particular case every particle carry two my_struct object
	 *
	 */

    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	// Here we define the boundary conditions of our problem
    size_t bc[2]={PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<2,float> g(0.01);

	// my_struct at position 0 in the aggregate
	constexpr int my_s1 = 0;

	// my_struct at position 1 in the aggregate
	constexpr int my_s2 = 1;

	//! \cond [vect create] \endcond

	vector_dist<2,float, aggregate<my_struct,my_struct>>
	vd(4096,domain,bc,g);

	std::cout << "HAS PACK: " << has_pack_agg<aggregate<my_struct,my_struct>>::result::value << std::endl;

	//! \cond [vect create] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop_ser Vector 4 property serialization
	 *
	 *
	 * ## Assign values to properties ##
	 *
	 * In this loop we assign position to particles and we fill the two my_struct
	 * that each particle contain. As demostration the first my_struct is filled
	 *  with the string representation of the particle coordinates. The second my struct
	 *   is filled with the string representation of the particle position multiplied by 2.0.
	 * The the vectors of the two my_struct are filled respectively with the sequence
	 * 1,2,3 and 1,2,3,4
	 *
	 *
	 *
	 * \snippet Vector/4_complex_prop/main_ser.cpp vect assign
	 *
	 */

	//! \cond [vect assign] \endcond

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		// we define x, assign a random position between 0.0 and 1.0
		vd.getPos(p)[0] = (float)rand() / RAND_MAX;

		// we define y, assign a random position between 0.0 and 1.0
		vd.getPos(p)[1] = (float)rand() / RAND_MAX;

		// Get the particle position as point
		Point<2,float> pt = vd.getPos(p);

		// create a C string from the particle coordinates
		// and copy into my struct
		vd.getProp<my_s1>(p).size = 32;
		vd.getProp<my_s1>(p).ptr = new char[32];
		strcpy(vd.getProp<my_s1>(p).ptr,pt.toString().c_str());

		// create a C++ string from the particle coordinates
		vd.getProp<my_s1>(p).str = std::string(pt.toString());

		vd.getProp<my_s1>(p).v.add(1);
		vd.getProp<my_s1>(p).v.add(2);
		vd.getProp<my_s1>(p).v.add(3);

		pt = pt * 2.0;

		// create a C string from the particle coordinates multiplied by 2.0
		// and copy into my struct
		vd.getProp<my_s2>(p).size = 32;
		vd.getProp<my_s2>(p).ptr = new char[32];
		strcpy(vd.getProp<my_s2>(p).ptr,pt.toString().c_str());

		// create a C++ string from the particle coordinates
		vd.getProp<my_s2>(p).str = std::string(pt.toString());

		vd.getProp<my_s2>(p).v.add(1);
		vd.getProp<my_s2>(p).v.add(2);
		vd.getProp<my_s2>(p).v.add(3);
		vd.getProp<my_s2>(p).v.add(4);

		// next particle
		++it;
	}

	//! \cond [vect assign] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop_ser Vector 4 property serialization
	 *
	 *
	 * ## Mapping and ghost_get ##
	 *
	 * Particles are redistributed across processors and we also synchronize the ghost
	 *
	 * \see \ref e0_s_map
	 *
	 * \see \ref e1_part_ghost
	 *
	 * \snippet Vector/4_complex_prop/main_ser.cpp vect map ghost
	 *
	 */

	//! \cond [vect map ghost] \endcond

	// Particles are redistribued across the processors
	vd.map();

	// Synchronize the ghost
	vd.ghost_get<my_s1,my_s2>();

	//! \cond [vect map ghost] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop_ser Vector 4 property serialization
	 *
	 *
	 * ## Output and VTK visualization ##
	 *
	 * Vector with complex properties can be still be visualized, because unknown properties are
	 * automatically excluded
	 *
	 * \see \ref e0_s_vis_vtk
	 *
	 * \snippet Vector/4_complex_prop/main.cpp vtk
	 *
	 */

	//! \cond [vtk] \endcond

	vd.write("particles");

	//! \cond [vtk] \endcond

	/*!
	 *
	 * \page Vector_4_complex_prop_ser Vector 4 property serialization
	 *
	 * ## Print 4 particles in the ghost area ##
	 *
	 * Here we print that the first 4 particles to show that the two my_struct contain the
	 * right information
	 *
	 * \snippet Vector/4_complex_prop/main_ser.cpp print ghost info
	 *
	 */

	//! \cond [print ghost info] \endcond

	size_t fg = vd.size_local();

	Vcluster & v_cl = create_vcluster();

	// Only the master processor print
	if (v_cl.getProcessUnitID() == 0)
	{
		// Print 4 particles
		for ( ; fg < vd.size_local()+4 ; fg++)
		{
			// Print my struct1 information
			std::cout << "my_struct1:" << std::endl;
			std::cout << "C-string: " << vd.getProp<my_s1>(fg).ptr << std::endl;
			std::cout << "Cpp-string: " << vd.getProp<my_s1>(fg).str << std::endl;

			for (size_t i = 0 ; i < vd.getProp<my_s1>(fg).v.size() ; i++)
				std::cout << "Element: " << i << "   " << vd.getProp<my_s1>(fg).v.get(i) << std::endl;

			// Print my struct 2 information
			std::cout << "my_struct2" << std::endl;
			std::cout << "C-string: " << vd.getProp<my_s2>(fg).ptr << std::endl;
			std::cout << "Cpp-string: " << vd.getProp<my_s2>(fg).str << std::endl;

			for (size_t i = 0 ; i < vd.getProp<my_s2>(fg).v.size() ; i++)
				std::cout << "Element: " << i << "   " << vd.getProp<my_s2>(fg).v.get(i) << std::endl;
		}
	}

	//! \cond [print ghost info] \endcond

	/*!
	 * \page Vector_4_complex_prop_ser Vector 4 property serialization
	 *
	 * ## Finalize ## {#finalize}
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet Vector/4_complex_prop/main_ser.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_4_complex_prop_ser Vector 4 property serialization
	 *
	 * # Full code # {#code}
	 *
	 * \include Vector/4_complex_prop/main_ser.cpp
	 *
	 */
}
