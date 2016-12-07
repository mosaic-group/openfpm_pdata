#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "VCluster.hpp"

/*
 *
 * \page VCluster_0_simple Using Vcluster to communicate across processors
 *
 *
 * [TOC]
 *
 * ## Simple example
 * 
 * This example show several basic functionalities of VCluster
 * 
 * 
 */
int main(int argc, char* argv[])
{
	/*
	 *
	 * \page VCluster_0_simple Using Vcluster to communicate across processors
	 *
	 *
	 * ## Initialization
	 *
	 * Before using any functionality the library must be initialized
	 *
	 * \snippet Vcluster/0_simple/main.cpp initialization
	 *
	 */

	//! \cond [initialization] \endcond

	openfpm_init(&argc,&argv);
	
	//! \cond [initialization] \endcond

	/*
	 *
	 * \page VCluster_0_simple Using Vcluster to communicate across processors
	 *
	 * ### Initialization of Vcluster
	 *
	 * Because in general our program is parallel we have more than one processors. With
	 * the function getProcessingUnits we can querry how many processors are involved in
	 * out computation
	 *
	 * \snippet Vcluster/0_simple/main.cpp initialization
	 *
	 */
	Vcluster & v_cl = create_vcluster();
	long int N_prc = v_cl.getProcessingUnits();

	/*
	 *
	 * \page VCluster_0_simple Using Vcluster to communicate across processors
	 *
	 *
	 * ### min max, sum
	 *
	 * with the function getProcessUnitID we can the the id of the processor executing
	 * the function. This function is equivalent to the MPI rank function.
	 * Vcluster provide several high and low level functionalities. One is max that
	 * return the maximum value across processors. There is also the function min
	 * and sum that return respectively the sum and the minimum across processors.
	 * All these operations are asynchronous, in order to get the result the function
	 * execute must be used. In our example the processor 0 print we also print the value
	 * but can be easily verified that also the other processors has the same value.
	 *
	 * \snippet Vcluster/0_simple/main.cpp max calc
	 *
	 */

	//! \cond [max calc] \endcond

	long int id = v_cl.getProcessUnitID();

	v_cl.max(id);
	v_cl.execute();
	if (v_cl.getProcessUnitID() == 0)
		std::cout << "Maximum processor rank: " << id << "\n";

	//! \cond [max calc] \endcond

	/*
	 *
	 * \page VCluster_0_simple Using Vcluster to communicate across processors
	 *
	 * We sum all the processor ranks the result should be that should
	 * be \$\frac{(n-1)n}{2}\$, only processor 0 print on terminal
	 *
	 *
	 * \snippet Vcluster/0_simple/main.cpp sum calc
	 *
	 */

	//! \cond [sum calc] \endcond

	size_t id2 = v_cl.getProcessUnitID();

	v_cl.sum(id2);
	v_cl.execute();
	if (v_cl.getProcessUnitID() == 0)
		std::cout << "Sum of all processors rank: " << id2 << "\n";

	//! \cond [sum calc] \endcond

	/*
	 *
	 * \page VCluster_0_simple Using Vcluster to communicate across processors
	 *
	 * We sum all the processor ranks the result should be that should
	 * be \$\frac{(n-1)n}{2}\$, only processor 0 print on terminal
	 *
	 *
	 * \snippet Vcluster/0_simple/main.cpp max calc
	 *
	 */

	long int id3 = v_cl.getProcessUnitID();
	openfpm::vector<long int> v;
	
	v_cl.allGather(id3,v);
	v_cl.execute();
	
	if (v_cl.getProcessUnitID() == 0)
	{
		std::cout << "Collected ids: ";
		for(size_t i = 0 ; i < v.size() ; i++)
			std::cout << " " << v.get(i) << " ";

		std::cout << "\n";
	}

	//
	// ### WIKI 5 ###
	//
	// we can also send messages to specific processors, with the condition that the receiving
	// processors know we want to communicate with them, if you are searching for a more
	// free way to communicate where the receiving processors does not know which one processor
	// want to communicate with us, see the example 1_dsde
	//

	std::stringstream ss_message_1;
	std::stringstream ss_message_2;
	ss_message_1 << "Hello from " << std::setw(8) << v_cl.getProcessUnitID() << "\n";
	ss_message_2 << "Hello from " << std::setw(8) << v_cl.getProcessUnitID() << "\n";
	std::string message_1 = ss_message_1.str();
	std::string message_2 = ss_message_2.str();
	size_t msg_size = message_1.size();
	
	// Processor 0 send to processors 1,2 , 1 to 2,1, 2 to 0,1

	v_cl.send(((id3+1)%N_prc + N_prc)%N_prc,0,message_1.c_str(),msg_size);
	v_cl.send(((id3+2)%N_prc + N_prc)%N_prc,0,message_2.c_str(),msg_size);

	openfpm::vector<char> v_one;
	v_one.resize(msg_size);
	openfpm::vector<char> v_two(msg_size);
	v_two.resize(msg_size);

	// Processor 0 receive from 1,2 ...

	v_cl.recv(((id3-1)%N_prc + N_prc)%N_prc,0,(void *)v_one.getPointer(),msg_size);
	v_cl.recv(((id3-2)%N_prc + N_prc)%N_prc,0,(void *)v_two.getPointer(),msg_size);
	v_cl.execute();

	if (v_cl.getProcessUnitID() == 0)
	{
		for (size_t i = 0 ; i < msg_size ; i++)
			std::cout << v_one.get(i);

		for (size_t i = 0 ; i < msg_size ; i++)
			std::cout << v_two.get(i);
	}

	//
	// ### WIKI 5 ###
	//
	// we can also do what we did before in one shot
	//

	id = v_cl.getProcessUnitID();
	id2 = v_cl.getProcessUnitID();
	id3 = v_cl.getProcessUnitID();
	v.clear();

	// convert the string into a vector

	openfpm::vector<char> message_1_v(msg_size);
	openfpm::vector<char> message_2_v(msg_size);

	for (size_t i = 0 ; i < msg_size ; i++)
		message_1_v.get(i) = message_1[i];

	for (size_t i = 0 ; i < msg_size ; i++)
		message_2_v.get(i) = message_2[i];

	v_cl.max(id);
	v_cl.sum(id2);
	v_cl.allGather(id3,v);

	// in the case of vector we have special functions that avoid to specify the size
	v_cl.send(((id+1)%N_prc + N_prc)%N_prc,0,message_1_v);
	v_cl.send(((id+2)%N_prc + N_prc)%N_prc,0,message_2_v);
	v_cl.recv(((id-1)%N_prc + N_prc)%N_prc,0,v_one);
	v_cl.recv(((id-2)%N_prc + N_prc)%N_prc,0,v_two);
	v_cl.execute();

	if (v_cl.getProcessUnitID() == 0)
	{
		std::cout << "Maximum processor rank: " << id << "\n";
		std::cout << "Sum of all processors rank: " << id << "\n";

		std::cout << "Collected ids: ";
		for(size_t i = 0 ; i < v.size() ; i++)
			std::cout << " " << v.get(i) << " ";

		std::cout << "\n";

		for (size_t i = 0 ; i < msg_size ; i++)
			std::cout << v_one.get(i);

		for (size_t i = 0 ; i < msg_size ; i++)
			std::cout << v_two.get(i);
	}

	openfpm_finalize();
}
