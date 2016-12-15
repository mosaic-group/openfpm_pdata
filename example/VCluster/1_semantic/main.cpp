#include "Grid/grid_dist_id.hpp"
#include "data_type/aggregate.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "VCluster/VCluster.hpp"

/*!
 *
 * \page VCluster_1_semantic Using Vcluster for Dynamic Sparse Data Exchange
 *
 * # Dynamic Sparse Data Exchange
 * 
 * Dynamic Sparse Data Exchange or DSDE, is a typical point to point communication in which
 * senders know to which processor to receive, but receivers has not knowledge about from
 * where they are receiving. OpenFPM use the NBX method or Non blocking consensus exchange.
 * (Said without bombastic world each processor wait for incoming messages. Pretty basic achivement
 * and technique in standard server programming, pictured bombastic and incredible discovery in MPI)
 * 
 */

#define N_NEXT 3
#define N_NUMBERS 5

int main(int argc, char* argv[])
{
	/*!
	 *
	 * \page VCluster_1_semantic using Vcluster for Dynamic Sparse Data Exchange
	 *
	 *
	 * ## Initialization
	 *
	 * Before using any functionality the library must be initialized. After initialization we can create
	 * the Vcluster object
	 *
	 * \snippet VCluster/0_simple/main.cpp initialization
	 *
	 */

	//! \cond [initialization] \endcond

	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();
	
	//! \cond [initialization] \endcond

	/*!
	 *
	 * \page VCluster_1_semantic Using Vcluster for Dynamic Sparse Data Exchange
	 *
	 * ## Dynamic Sparse Data Exchange
	 *
	 * To do dynamic sparse data exchange, each processor fill a send processor list
	 * and create a message for each processor. In this case the message will be a complex
	 * object. OpenFPM use the capability to serialize complex object into sequence of byte
	 * to send over the network and de-serializa or re-assemble the object into another
	 * processors. In this case the complex object is a list of double numbers. At the end
	 * of the example each processor print what it received
	 *
	 * \snippet VCluster/1_semantic/main.cpp ssendrecv
	 *
	 */

	//! \cond [ssendrecv] \endcond

	// id of the processor calling this function
	long int proc_id = v_cl.getProcessUnitID();

	// number of processors executing this program
	long int n_proc = v_cl.getProcessingUnits();

	// List of processors we communicate with
	openfpm::vector<size_t> prc_send;

	// For each processor we want to send a vector of doubles
	// in this case each processor send N_NEXT vectors.
	// In general we can think to openfpm::vector<T> as a set of objects T.
	// where we want to send one object T to each processor in out sending list.
	// In this case T is a list of double or (openfpm::vector<double>)
	openfpm::vector<openfpm::vector<double>> messages_send(N_NEXT);

	// Here we prepare the senbding buffer
	for (size_t i = 0, m = 0 ; i < N_NEXT ; i++, m++)
	{
		// create the sending processor list
		prc_send.add(openfpm::math::positive_modulo(proc_id + i + 1,n_proc));

		// Fill with soma data the vectors
		for (size_t j = 0 ; j < N_NUMBERS ; j++)
			messages_send.get(m).add(j+N_NUMBERS*proc_id);

	}

	// Buffer that receive messages
	openfpm::vector<double> messages_recv2;

	// List of processor from which we receove
	openfpm::vector<size_t> prc_recv2;

	// number of elements we receive from each processors
	openfpm::vector<size_t> sz_recv2;

	v_cl.SSendRecv(messages_send,messages_recv2,prc_send,prc_recv2,sz_recv2);

	// here each processor print the received message

	std::cout << "Processor " << proc_id << " received ";
	for (size_t i = 0 ; i < messages_recv2.size() ; i++)
		std::cout << messages_recv2.get(i) << "   ";

	std::cout << std::endl;

	//! \cond [ssendrecv] \endcond

	/*!
	 * \page VCluster_1_semantic Using Vcluster to communicate across processors
	 *
	 * ## Finalize ##
	 *
	 *  At the very end of the program we have always to de-initialize the library
	 *
	 * \snippet VCluster/1_semantic/main.cpp finalize
	 *
	 */

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond
}
