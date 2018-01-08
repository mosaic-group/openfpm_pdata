/*!
 * \page PS-CMA-ES Particle swarm CMA Evolution strategy
 *
 *
 * [TOC]
 *
 *
 * # Optimization # {#Opti_cma_es}
 *
 *
 * In this example we show how to code PS-CMA-ES. This is just a simple variation to the
 * CMA-ES, where you have multiple CMA-ES running. The the best solution across them is
 * used to produce a drift velocity toward that point.
 *
 *
 *
 */

#include "Vector/vector_dist.hpp"
#include "Eigen/Dense"
#include <Eigen/Eigenvalues>
#include <Eigen/Jacobi>
#include <limits>
#include "Vector/vector_dist.hpp"

constexpr int dim = 2;
// set this to 4+std::floor(3*log(dim))
constexpr int lambda = 6;
constexpr int mu = lambda/2;

constexpr int sigma = 0;
constexpr int sample = 1;
constexpr int Cov_m = 2;
constexpr int B = 3;
constexpr int D = 4;
constexpr int Zeta = 5;
constexpr int path_s = 6;
constexpr int path_c = 7;
constexpr int ord = 8;

const double c_m = 1.0;

double mu_eff = 1.0;
double cs = 1.0;
double cc = 1.0;
double c1 = 1.0;
double c_mu = 0.5;
double chiN;
double d_amps = 1.0;
double stop_fitness = 1.0;
int eigeneval = 0;

typedef vector_dist<dim,double, aggregate<double,
										 Eigen::VectorXd[lambda],
										 Eigen::MatrixXd,
										 Eigen::MatrixXd,
										 Eigen::DiagonalMatrix<double,Eigen::Dynamic>,
										 Eigen::VectorXd[lambda],
										 Eigen::VectorXd,
										 Eigen::VectorXd,
										 int[lambda]> > particle_type;

double f(Eigen::VectorXd & v)
{
	double ret = 0.0;

	ret = v.transpose()*v;

	return ret;
}

double generateGaussianNoise(double mu, double sigma)
{
	static const double epsilon = std::numeric_limits<double>::min();
	static const double two_pi = 2.0*3.14159265358979323846;

	thread_local double z1;
	thread_local double generate;
	generate = !generate;

	if (!generate)
	{return z1 * sigma + mu;}

	double u1, u2;
	do
	{
	   u1 = rand() * (1.0 / RAND_MAX);
	   u2 = rand() * (1.0 / RAND_MAX);
	}
	while ( u1 <= epsilon );

	double z0;
	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}

template<unsigned int dim>
Eigen::VectorXd generateGaussianVector()
{
	Eigen::VectorXd tmp;
	tmp.resize(dim);

	for (size_t i = 0 ; i < dim ; i++)
	{
		tmp(i) = generateGaussianNoise(0,1);
	}

	return tmp;
}

template<unsigned int dim>
void fill_vector(double (& f)[dim], Eigen::VectorXd & ev)
{
	for (size_t i = 0 ; i < dim ; i++)
	{ev(i) = f[i];}
}

void fill_vector(const double * f, Eigen::VectorXd & ev)
{
	for (size_t i = 0 ; i < ev.size() ; i++)
	{ev(i) = f[i];}
}

struct fun_index
{
	double f;
	int id;

	bool operator<(const fun_index & tmp)
	{
		return f < tmp.f;
	}
};

double wm[mu];

void init_weight()
{
	for (size_t i = 0 ; i < mu ; i++)
	{wm[i] = log(mu+0.5) - log(i+1);}

	double tot = 0.0;

	for (size_t i = 0 ; i < mu ; i++)
	{tot += wm[i];}

	double sum = 0.0;
	double sum2 = 0.0;

	for (size_t i = 0 ; i < mu ; i++)
	{
		wm[i] /= tot;
		sum += wm[i];
		sum2 += wm[i]*wm[i];
	}

	// also set mu_eff
    mu_eff=sum*sum/sum2;

}

double weight(int i)
{
	return wm[i];
}

void cma_step(particle_type & vd, int step,  double & best, int & best_i)
{
	best = std::numeric_limits<double>::max();
	auto it = vd.getDomainIterator();

	Eigen::VectorXd mean_x_new(dim);
	Eigen::VectorXd mean_x_old(dim);
	Eigen::VectorXd mean_z(dim);

	openfpm::vector<fun_index> f_obj(lambda);

	int counteval = step*lambda;

	while (it.isNext())
	{
		auto p = it.get();

		// fill the mean vector;

		fill_vector(vd.getPos(p),mean_x_old);

		for (size_t j = 0 ; j < lambda ; j++)
		{
			vd.getProp<Zeta>(p)[j] = generateGaussianVector<dim>();
			vd.getProp<sample>(p)[j] = mean_x_old + vd.getProp<sigma>(p)*vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<Zeta>(p)[j];
			f_obj.get(j).f = f(vd.getProp<sample>(p)[j]);
			f_obj.get(j).id = j;
		}

		f_obj.sort();

		for (size_t j = 0 ; j < lambda ; j++)
		{vd.getProp<ord>(p)[j] = f_obj.get(j).id;}

		// Calculate weighted mean

		mean_x_new.setZero();
		mean_z.setZero();
		for (size_t j = 0 ; j < mu ; j++)
		{
			mean_x_new += weight(j)*vd.getProp<sample>(p)[vd.getProp<ord>(p)[j]];
			mean_z += weight(j)*vd.getProp<Zeta>(p)[vd.getProp<ord>(p)[j]];
		}

		vd.getProp<path_s>(p) = vd.getProp<path_s>(p)*(1.0 - cs) + sqrt(cs*(2.0-cs)*mu_eff)*vd.getProp<B>(p)*mean_z;

		double hsig = vd.getProp<path_s>(p).norm()/(1-pow(1-cs,2*counteval/lambda))/dim < 2.0 + 4.0/(dim+1);

		vd.getProp<path_c>(p) = (1-cc)*vd.getProp<path_c>(p) + hsig * sqrt(cc*(2-cc)*mu_eff)*(vd.getProp<B>(p)*vd.getProp<D>(p)*mean_z);

		// Adapt covariance matrix C
		vd.getProp<Cov_m>(p) = (1-c1-c_mu)*vd.getProp<Cov_m>(p) +
				            c1*(vd.getProp<path_c>(p)*vd.getProp<path_c>(p).transpose() + (1-hsig)*cc*(2-cc)*vd.getProp<Cov_m>(p));

		for (size_t i = 0 ; i < mu ; i++)
		{vd.getProp<Cov_m>(p) += c_mu*(vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<Zeta>(p)[vd.getProp<ord>(p)[i]])*weight(i)*
			                          (vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<Zeta>(p)[vd.getProp<ord>(p)[i]]).transpose();
		}

		//Adapt step-size sigma
		vd.getProp<sigma>(p) = vd.getProp<sigma>(p)*exp((cs/d_amps)*(vd.getProp<path_s>(p).norm()/chiN - 1));

		std::cout << vd.getProp<sigma>(p) << std::endl;

		// Update B and D from C

		if (counteval - eigeneval > lambda/(c1+c_mu)/dim/10)
		{
			eigeneval = counteval;
			//vd.getProp<Cov_m>(p) = (vd.getProp<Cov_m>(p)+vd.getProp<Cov_m>(p).transpose()) / 2.0; // enforce symmetry

			// Eigen decomposition
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver;
			eig_solver.compute(vd.getProp<Cov_m>(p));

			for (size_t i = 0 ; i < eig_solver.eigenvalues().size() ; i++)
			{vd.getProp<D>(p).diagonal()[i] = sqrt(eig_solver.eigenvalues()[i]);}
			vd.getProp<B>(p) = eig_solver.eigenvectors();
		}

	    // Break, if fitness is good enough or condition exceeds 1e14, better termination methods are advisable
//	    if (f_obj.get(0).f <= stop_fitness || vd.getProp<D>(p).diagonal().maxCoeff() > 1e7 * vd.getProp<D>(p).diagonal().minCoeff())
//	    {break;}

	    // Escape flat fitness, or better terminate?

	    if (f_obj.get(0).f == f_obj.get(std::ceil(0.7*lambda)).f )
	    {
	    	vd.getProp<sigma>(p) = vd.getProp<sigma>(p)*exp(0.2+cs/d_amps);
	    	std::cout << "warning: flat fitness, consider reformulating the objective";
	    }

	    // Copy the new mean as position of the particle
	    for (size_t i = 0 ; i < dim ; i++)
	    {vd.getPos(p)[i] = mean_x_new(i);}

	    if (best > f_obj.get(0).f)
	    {
	    	best = f_obj.get(0).f;
	    	best_i = p.getKey();
	    }
	    std::cout << "Best solution: " << f_obj.get(0).f << "   " << vd.getProp<sigma>(p) << std::endl;

		++it;
	}
}

void broadcast_best_solution(particle_type & vd, openfpm::vector<double> & best_sol, double best, size_t best_i)
{
	best_sol.resize(dim);
	auto & v_cl = create_vcluster();

	double best_old = best;
	v_cl.min(best);
	v_cl.execute();

	size_t rank;
	if (best_old == best)
	{
		rank = v_cl.getProcessUnitID();

		// we own the minimum and we decide who broad cast
		v_cl.min(rank);
		v_cl.execute();

		if (rank == v_cl.getProcessUnitID())
		{
			for (size_t i = 0 ; i < dim ; i++)
			{best_sol.get(i) = vd.getPos(best_i)[i];}
		}
	}
	else
	{
		rank = std::numeric_limits<size_t>::max();

		// we do not own  decide who broad cast
		v_cl.min(rank);
		v_cl.execute();
	}

	// now we broad cast the best solution across processors

	v_cl.Bcast(best_sol,rank);
	v_cl.execute();
}

void write_million_point_on_file()
{
	std::ofstream output("rnd_output");

	for (size_t i = 0 ; i < 1000000 ; i++)
	{
		double rnd = generateGaussianNoise(0.0,1.0);
		output << std::setprecision(15) << rnd << std::endl;
	}

	output.close();
}


void create_rotmat(Eigen::VectorXd & S,Eigen::VectorXd & T, Eigen::MatrixXd & R)
{
	Eigen::VectorXd S_work(dim);
	Eigen::VectorXd T_work(dim);
	Eigen::VectorXd S_sup(dim);
	Eigen::VectorXd T_sup(dim);

	Eigen::MatrixXd R_tar(dim,dim);
	Eigen::MatrixXd R_tmp(dim,dim);
	Eigen::MatrixXd R_sup(dim,dim);
	double G_S,G_C;
	Eigen::MatrixXd S_tmp(2,2);
	Eigen::MatrixXd T_tmp(2,2);
	int p,q,i;

	S_work = S;
	T_work = T;

	R.setIdentity();
	R_tar = R;
	R_tmp = R;

	for (p = dim - 2; p >= 0 ; p -= 1)
	{

		for (q = dim - 1 ; q >= p+1 ; q-= 1)
		{
			T_tmp(0) = T_work(p);
			T_tmp(1) = T_work(q);
			S_tmp(0) = S_work(p);
			S_tmp(1) = S_work(q);

			// Perform Givens Rotation on start vector

			Eigen::JacobiRotation<double> G;
			G.makeGivens(S_tmp(0), S_tmp(1));

			// Check direction of rotation
			double sign = 1.0;
			if (S_tmp(1) > 0.0)
			{sign = -1.0;}

			// Build a Rotation Matrix out of G_C and G_S
			R_tmp.setIdentity();
			R_tmp(p,p) = sign*G.c();
			R_tmp(q,q) = sign*G.c();
			R_tmp(p,q) = sign*G.s();
			R_tmp(q,p) = sign*-G.s();

			// Rotate start vector and update R
			// S_work = R_tmp*S_work

			S_work = R_tmp*S_work;
			// R = R_tmp*R
			R = R_tmp*R;

			// Perform Givens Rotation on target vector

			G.makeGivens(T_tmp(0), T_tmp(1));

			sign = 1.0;
			if (T_tmp(1) < 0.0)
			{sign = -1.0;}

			R_tmp.setIdentity();
			R_tmp(p,p) = sign*G.c();
			R_tmp(q,q) = sign*G.c();
			R_tmp(p,q) = sign*G.s();
			R_tmp(q,p) = sign*-G.s();

			// Rotate target vector and update R_tar

			T_work = R_tmp*T_work;
			R_tar = R_tmp*R_tar;
		}
	}

	R = R_tar.transpose()*R;

	// Check the rotation

	Eigen::VectorXd Check(dim);
	Check = R*S;
}

void rotate_covariant_matrix_and_bias(particle_type & vd, const openfpm::vector<double> & best_sol)
{
	auto it = vd.getDomainIterator();

	Eigen::VectorXd eigen_v(dim);
	Eigen::VectorXd best_sol_v(dim);
	fill_vector(&best_sol.get(0),best_sol_v);
	Eigen::VectorXd pos;

	// Calculate the target
	Eigen::MatrixXd R(dim,dim);
	// Target
	Eigen::VectorXd T(dim);

	while (it.isNext())
	{
		auto p = it.get();

		fill_vector(vd.getPos(p),pos);

		int max_i;
		vd.getProp<D>(p).diagonal().maxCoeff(&max_i,&max_i);

		eigen_v = vd.getProp<B>(p).col(max_i);

		////////// DEBUG ///////////////

		Eigen::VectorXd debug(dim);
		Eigen::VectorXd debug2(dim);

		debug2 = vd.getProp<Cov_m>(p)*debug;

		///////////////////////////////////
		T = best_sol_v - pos;

		// Now we rotate
		create_rotmat(eigen_v,T,R);

		++it;
	}
}


int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

//	write_million_point_on_file();
//	return 0;

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<dim,double> domain;

	for (size_t i = 0 ; i < dim ; i++)
	{
		domain.setLow(i,0.0);
		domain.setHigh(i,1.0);
	}

	// Here we define the boundary conditions of our problem
	size_t bc[dim];
	for (size_t i = 0 ; i < dim ; i++)
    {bc[i] = NON_PERIODIC;};

	// extended boundary around the domain, and the processor domain
	Ghost<dim,double> g(0.0);

    particle_type vd(1,domain,bc,g);

    // Initialize constants

    stop_fitness = 1e-10;
    size_t stopeval = 1e3*dim*dim;

    // Strategy parameter setting: Selection
    init_weight();

    // Strategy parameter setting: Adaptation
    cc = (4.0 + mu_eff/dim) / (dim+4.0 + 2.0*mu_eff/dim);
    cs = (mu_eff+2) / (dim+mu_eff+5);
    c1 = 2.0 / ((dim+1.3)*(dim+1.3)+mu_eff);
    c_mu = std::min(1-c1, 2 * (mu_eff-2+1/mu_eff) / ((dim+2)*(dim+2)+mu_eff));
    d_amps = 1 + 2*std::max(0.0, sqrt((mu_eff-1)/(dim+1))-1) + cs;

    chiN = sqrt(dim)*(1.0-1.0/(4.0*dim)+1.0/(21.0*dim*dim));

	//! \cond [assign position] \endcond

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		for (size_t i = 0 ; i < dim ; i++)
		{
			// we define x, assign a random position between 0.0 and 1.0
			vd.getPos(p)[i] = 1.0/*(double)rand() / RAND_MAX*/;
		}

		vd.getProp<sigma>(p) = 0.5;

		// Initialize the covariant Matrix,B and D to identity

		vd.getProp<D>(p).resize(dim);
		for (size_t i = 0 ; i < vd.getProp<D>(p).diagonal().size() ; i++)
		{vd.getProp<D>(p).diagonal()[i] = 1.0;}
		vd.getProp<B>(p).resize(dim,dim);
		vd.getProp<B>(p).setIdentity();
		vd.getProp<Cov_m>(p) = vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<D>(p)*vd.getProp<B>(p);
		vd.getProp<path_s>(p).resize(dim);
		vd.getProp<path_s>(p).setZero(dim);
		vd.getProp<path_c>(p).resize(dim);
		vd.getProp<path_c>(p).setZero(dim);

		// next particle
		++it;
	}

	double best = 0.0;
	int best_i = 0;

	openfpm::vector<double> best_sol;
	// now do several iteration

	for (size_t i = 0 ; i < 100 ; i++)
	{
		// sample offspring
		cma_step(vd,i+1,best,best_i);

		// Find the best point across processors
		broadcast_best_solution(vd,best_sol,best,best_i);

		// Generate the rotational Matrix
		rotate_covariant_matrix_and_bias(vd,best_sol);

		//create_rotmat();

		// Rotate the Covariant MAtrix
	}

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_0_simple Vector 0 simple
	 *
	 * ## Full code ## {#code_e0_sim}
	 *
	 * \include Vector/0_simple/main.cpp
	 *
	 */
}
