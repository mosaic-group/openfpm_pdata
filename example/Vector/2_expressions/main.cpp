#include "Vector/vector_dist.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

/*
 * ### WIKI 1 ###
 *
 * ## Simple example to show how expression work with vector
 * 
 * This example show how to use expressions to do complex calculation on particles
 * 
 * ### WIKI END ###
 * 
 */

constexpr int A = 0;
constexpr int B = 1;
constexpr int C = 2;
constexpr int D = 3;


// First we define the kernel function (an exponential in this case)
struct exp_kernel
{
	// variance of the exponential
	float sigma;

	// Constructor require to define the variance
	exp_kernel(float sigma)
	:sigma(sigma)
	{}

	// The kernel function accept the position p, the position q
	inline Point<3,double> value(const Point<3,float> & p, const Point<3,float> & q,Point<3,double> & Pp,Point<3,double> & Pq)
	{
		// calculate the distance between p and q
		float dist = norm(p-q);

			
		return Pq * exp(- dist * dist / sigma);
	}
};


int main(int argc, char* argv[])
{
	//
	// ### WIKI 2 ###
	//
	// Here we Initialize the library, than we create a uniform random generator between 0 and 1 to to generate particles
	// randomly in the domain, we create a Box that define our domain, boundary conditions, and ghost
	//
	
	// initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with intervals from 0 to 1.0
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Here we define the boundary conditions of our problem
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<3,float> g(0.01);

	// Delta t
	double dt = 0.01;

	// distributed vectors
	vector_dist<3,double, aggregate<double,double,Point<3,double>,Point<3,double>> > vd(4096,domain,bc,g);
	vector_dist<3,double, aggregate<double> > vd2(4096,domain,bc,g);

	// Assign random position to the vector dist
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		vd.getPos(p)[0]  = (double)rand() / RAND_MAX;
		vd.getPos(p)[1]  = (double)rand() / RAND_MAX;
		vd.getPos(p)[2]  = (double)rand() / RAND_MAX;

		vd2.getPos(p)[0] = (double)rand() / RAND_MAX;
		vd2.getPos(p)[1]  = (double)rand() / RAND_MAX;
		vd2.getPos(p)[2]  = (double)rand() / RAND_MAX;

		++it;
	}

	// Redistribute the particles and synchronize the ghost
	vd.map();
	vd.map();

	//
	// ### WIKI ###
	//
	// Here we create alias for the properties of the vectors
	//

	// vA is an alias for the property A of the vector vd
	// vB is an alias for the property V of the vector vd
	// vC is ...
	auto vA = getV<A>(vd);
	auto vB = getV<B>(vd);
	auto vC = getV<C>(vd);
	auto vD = getV<D>(vd);
	auto vPOS = getV<PROP_POS>(vd);

	// same concept for the vector v2
	auto v2A = getV<0>(vd);
	auto v2POS = getV<PROP_POS>(vd);

	// Assign 1 to the property A
	vA = 1;

	// Equal the property A to B
	vB = vA;

	// Assign to the property C and for each component 4.0
	vC = 4;

	// Assign te position of the particle to the property D
	vD = vPOS;

	//
	// ### WIKI ###
	//
	// All these expression are applied point wise for each particles
	//
	// P.S.
	//
	// the * is a scalar product in case of vectors
	// simple multiplication in case of scalar
	//
	//

	// In vA we set the distance from the origin for each particle
	vA = norm(vPOS);

	//
	// For each particle p calculate the expression under
	//
	// NOTE sin(2.0 * vD)  and exp(5.0 * vD) are both vector
	//
	//              so this is a scalar product
	//                      |
	//                      V
	vA = vA + sin(2.0 * vD) * exp(1.2 * vD);
    //       |_____________________________|
    //                      |
    //                      V
    //              and this is a number
    //
    //


	// pmul indicate component-wise multiplication the same as .* in Matlab
	vC = pmul(vD,vD) * dt;

	// Normalization of the vector vD
	vD = vD / sqrt( vD * vD );

	// we write the result on file
	vd.write("output");

	// We can also apply the expression in a second moment like this
	Point<3,double> p0({0.5,0.5,0.5});

	// we cannot use p0 directly we have to create an expression from it
	auto p0_e = getVExpr(p0);

	// As soon as the second vector have the some number of particle as v we can apply the expr1 to the v2 with v2D
	auto expr1 = 1.0/2.0/sqrt(M_PI)*exp(-(v2POS-p0_e)*(v2POS-p0_e)/2.0);

	// here is executed on all particles of v2 and the result assigned to the property A of v2
	v2A = expr1;

	//
	// ### WIKI ###
	//
	// Each operator= produce a iteration across particles. It is possible to use the function assign
	// to execute multiple instructions in one cycle
	//
	// assign support an variable of istructions
	//

	assign(vA,1.0,               // vA = 1
		   vB,expr1,          // vB = 1.0/2.0/sqrt(M_PI)*exp(-(v2POS-p0_e)*(v2POS-p0_e)/2.0)
		   vC,vD/norm(vD),    // vC = vD/norm(vD)
		   vD,2.0*vC);        // vD = 2.0*vC

	// As soon as the vector v have the some number of particle as v2 we can apply the expr1 to v
	vA = expr1;

	//
	// ### WIKI ###
	//
	// we shown the simple usage of simple point-wise function like sin
	// and exp. Here we show the function applyKernel than encapsulate
	// the operation
	//
	// \sum_{q = Neighborhood p} Kernel(x_p,x_q,Up,Uq) * expression
	//
	// where x_p is the position of the particle p
	//       x_q is the position of the particle q
	//       U_p is the value of the property carried by the particle p
	//       U_q is the value of the property carried by the particle q
	//
	//

	exp_kernel ker(0.5);
	auto cl = vd.getCellList(0.1);

	// We are going to do some calculation that require ghost, so sync it
        vd.ghost_get<3>();

	// Here we create the

	// A renormalizaion of the vector vD is calculated on fly
	//
	// and than For each particle p
	//
	// \sum_{q = Neighborhood p} Ker(x_p,x_q,D_p,D_q) * \frac{D}{norm(D)}
	//
	//
	// P.S.
	//
	//  Cerefull here ApplyKernel is a special function if it act on property D we cannot write on property D 
	//  ( vD = ...  is an error and produce unpredictable result )
	//
	vC = applyKernel_in( vD / norm(vD) ,vd,cl,ker) + vD;

	// ok, because does not use applyKernel
	vD = vD / norm(vD);

	//
	// ### WIKI 9 ###
	//
	// Output the particle position for each processor
	//

	vd.deleteGhost();
	vd.write("output2");

	//
	// ### WIKI 10 ###
	//
	// Deinitialize the library
	//
	openfpm_finalize();
}
