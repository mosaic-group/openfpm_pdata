/*!
 *
 * \page PS_CMA_ES Particle swarm CMA-ES Evolution strategy
 *
 *
 * [TOC]
 *
 * # Optimization {#Opti_cma_es}
 *
 *
 * In this example we show how to code PS-CMA-ES. This is just a simple variation to the
 * CMA-ES, where you have multiple CMA-ES running. The the best solution across them is
 * used to produce a drift velocity toward that point.
 *
 * ## Introduction {#ps_cme_es}
 *
 * In this example we try to find the global optimum of a function. In particular we are
 * using the function F15 from the CEC 2005 benchmark test, to validate that PS-CMA-ES work.
 * This example contain multiple files:
 *
 * * f15_cec_const.hpp definitions of constants for the F15 function
 * * f15_cec_fun.hpp the function itself
 *
 * The function is quite complicated and for reference please refere to the function
 *  F15 "Hybrid Composition" in the CEC 2005 test. The function can be called with
 *  hybrid_composition<dim>(x) where dim is the dimensionality and x is the point
 *   where is evaluated the function. The dimensionality can go from 1 to 50.
 *
 * Considering to have a function \f$ f \f$ from \f$ \mathbb{R}^{dim} \f$ to \f$ \mathbb{R} \f$,
 * the algorithm use a set of particles to find in parallel the global optimum of a function.
 *  The algorithm rather than try to find the global optimum
 * sampling point randomly in the space, it uses a set of particles each of them having a gaussian
 * sampling distribution \f$ e^{\sigma \cdot x^tCx} \f$ with C a \f$ dim \cdot dim \f$ matrix.
 * At each step for each particle p **lambda** points are sampled using the sampling
 * distribution centered on the particle position. The covariant matrix and sigma are is subsequently
 * adjusted to favor sampling around the best sampled points. In order to do this the algorithm
 * need the eigen-value decomposition of \f$ C = B^{t}DB \f$ where \f$ D \f$ is a diagonal
 * Matrix and \f$ B \f$ is the Matrix of the Eigen-vector. In order to reduce or increase
 * the sampling area the sigma is instead used. The algorithm use the vector **path_s** to
 * detect stagnation of the particle movement, and use **path_c**(a transfomed version of **path_s**)
 * to refine the sampling covariant matrix from the fact that the particle is "moving" toward that
 * direction. PS-CMA-ES is just a variation in which every **N_pos** CMA-Es steps the CMA-ES is
 * sampling distribution and position is biased toward the best founded point across all independent
 * CMA-ES.
 *
 * Explain the CMA-ES algorithm in detail is out of the purpose of this tutorial example.
 * We will briefly go across the main step. For a full reference of the CMA-ES
 * algoritm please refers to <a href="https://arxiv.org/abs/1604.00772">this paper</a>.
 * While for PS-CMA-ES refers to <a href="http://mosaic.mpi-cbg.de/docs/Mueller2009a.pdf">this paper</a>.
 *
 *
 * ## Inclusions {#inclusions_and_constants}
 *
 * In this example we use a set of particles so we will use **vector_dist**, we will use
 * Eigen dense matrix. Because the standard dense matrix are not compatible with
 * the vector_dist we will use **EMatrix** that are simple wrapper to Eigen::Matrix
 * but compatible with vector_dist. Because EMatrix are compatible with all the
 * Eigen value functions we can use them in all Eigen functions. For CMA-ES algorithm
 * we also need Eigen-value eigen-vector decomposition and Jacobi-Rotation for the
 * particle-swarm part.
 *
 * \snippet Numerics/PS-CMA-ES/main.cpp ps_cma_es_inclusion
 *
 * PS-CMA-ES require several properties to be stored on the particles, some has been already
 * explained. Here we explain the others.
 *
 * * **Zeta** contain the lambda sampled points (before apply the covariant matrix tranformation)
 *        so it contain points samples on a gaussian of sigma=1 centered in zero
 *
 * * **ord** Contain the sequrnce if we want the lambda generated points in order from the best to
 *       the worst
 *
 * * **stop** If a flag that indicate that the CMA-ES reached some stop criteria
 *
 * * **fithist** It contain historical information about the particles to penalize them in case the
 *       go out of boundary. It is 1:1 taken from cmaes.m (production version)
 *       <a href="https://www.lri.fr/~hansen/cmaes_inmatlab.html">this paper</a> (or Google it)
 *
 * * **weight** Same concept of fithist other information to penalize particles going out of
 *       the boundary
 *
 * * **validfit** Same concept of fithist other information to penalize particles going out of
 *       the boundary
 *
 * * **xold** It contain the previous position of the particles used in several calculations
 *
 * * **last_restart** CMA-ES The CMA-ES sigma become very small the CMA-ES converged. At this point
 *       we can do two things, one is to stop the CMA-ES, the other is to restart-it to explore
 *       better the space. In case it restart. this parameter indicate at which iteration happen the
 *       last restart.
 *
 * * **iniphase** Same concept of fithist other information to penalize particles going out of
 *       the boundary
 *
 * * **xmean_st** This contain the new position of the particle it will be stored as particle position
 *       at the end of the CMA-ES step
 *
 * * **xmean_st** This contain the new position of the particle in a space where we do not apply the
 *       covariant transformation. (In practice is a weighted sum of the Zeta samples points)
 *
 * \snippet Numerics/PS-CMA-ES/main.cpp def_part_set
 *
 * ## Parameters {#ps_cma_par}
 *
 * CMA-ES and further PS-CMA-ES require some parameters in order to work. refers to the
 * papers above to have a full explanation, but here is a short one
 *
 * * **dim** Dimensionality of the test function
 *
 * * **lambda** number of sample points taken at each iteration by CMA-ES
 *              suggested to use \f$ 4+floor(3*log(dim)) \f$
 *
 * * **mu** only mu best points are considered to adapt the Covariant matrix
 *
 * * **psoWeight** How much the pso step bias the particle positions
 *
 * * **N_pso** Number of CMA-ES step before do a PSO step (200 give the possibility
 *   to the CMA-ES to explore the neighborhood and narrow down at least a funnel)
 *
 * * **stopTolX** stop criteria for CMA-ES. When the the sampling area is small enough
 *   stop
 *
 * * **StopToUpX** stop criteria is the sampling area become too big
 *
 * * **restart_cma** If the CMA-ES reach a stop criteria reinitialize and restart
 *
 * * **hist_size** size of the array fit_hist (default should be mainly fine)
 *
 */

//! [ps_cma_es_inclusion]

#define EIGEN_USE_LAPACKE
#include "Vector/vector_dist.hpp"
#include "DMatrix/EMatrix.hpp"
#include <Eigen/Eigenvalues>
#include <Eigen/Jacobi>
#include <limits>
#include "Vector/vector_dist.hpp"
#include <f15_cec_fun.hpp>
#include <boost/math/special_functions/sign.hpp>

//! [ps_cma_es_inclusion]

//! [parameters]

// PARAMETERS
constexpr int dim = 10;
// when you set dim set also lambda to to 4+std::floor(3*log(dim))
constexpr int lambda = 7;
constexpr int mu = lambda/2;
double psoWeight = 0.7;
// number of cma-step before pso step
int N_pso = 200;
double stopTolX = 2e-11;
double stopTolUpX = 2000.0;
int restart_cma = 1;
size_t max_fun_eval = 30000000;
constexpr int hist_size = 21;

// Convenient global variables (Their value is set after)
double mu_eff = 1.0;
double cs = 1.0;
double cc = 1.0;
double ccov = 1.0;
double chiN;
double d_amps = 1.0;
double stop_fitness = 1.0;
int eigeneval = 0;
double t_c = 0.1;
double b = 0.1;

//! [parameters]

//! [def_part_set]

//////////// definitions of the particle set

constexpr int sigma = 0;
constexpr int Cov_m = 1;
constexpr int B = 2;
constexpr int D = 3;
constexpr int Zeta = 4;
constexpr int path_s = 5;
constexpr int path_c = 6;
constexpr int ord = 7;
constexpr int stop = 8;
constexpr int fithist = 9;
constexpr int weight = 10;
constexpr int validfit = 11;
constexpr int xold = 12;
constexpr int last_restart = 13;
constexpr int iniphase = 14;
constexpr int xmean_st = 15;
constexpr int meanz_st = 16;

typedef vector_dist<dim,double, aggregate<double,
										 Eigen::MatrixXd,
										 Eigen::MatrixXd,
										 Eigen::DiagonalMatrix<double,Eigen::Dynamic>,
										 Eigen::VectorXd[lambda],
										 Eigen::VectorXd,
										 Eigen::VectorXd,
										 int[lambda],
										 int,
										 double [hist_size],
										 double [dim],
										 double,
										 Eigen::VectorXd,
										 int,
										 bool,
										 Eigen::VectorXd,
										 Eigen::VectorXd> > particle_type;

//! [def_part_set]


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
	z0 = sqrt(-2.0 * log(u2)) * cos(two_pi * u1);
	z1 = sqrt(-2.0 * log(u2)) * sin(two_pi * u1);
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
	{wm[i] = log(double(mu)+1.0) - log(double(i)+1.0);}

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

double weight_sample(int i)
{
	return wm[i];
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
			double z;
			G.makeGivens(S_tmp(0), S_tmp(1),&z);

			// Check direction of rotation
			double sign = 1.0;
			if (z < 0.0)
			{sign = -1.0;}

			// Build a Rotation Matrix out of G_C and G_S
			R_tmp.setIdentity();
			R_tmp(p,p) = sign*G.c();
			R_tmp(q,q) = sign*G.c();
			R_tmp(p,q) = sign*-G.s();
			R_tmp(q,p) = sign*G.s();

			// Rotate start vector and update R
			// S_work = R_tmp*S_work

			S_work = R_tmp*S_work;
			// R = R_tmp*R
			R = R_tmp*R;

			// Perform Givens Rotation on target vector

			G.makeGivens(T_tmp(0), T_tmp(1),&z);

			sign = 1.0;
			if (z < 0.0)
			{sign = -1.0;}

			R_tmp.setIdentity();
			R_tmp(p,p) = sign*G.c();
			R_tmp(q,q) = sign*G.c();
			R_tmp(p,q) = sign*-G.s();
			R_tmp(q,p) = sign*G.s();

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

void updatePso(openfpm::vector<double> & best_sol,
			   double sigma,
			   Eigen::VectorXd & xmean,
			   Eigen::VectorXd & xold,
			   Eigen::MatrixXd & B,
			   Eigen::DiagonalMatrix<double,Eigen::Dynamic> & D,
			   Eigen::MatrixXd & C_pso)
{
	Eigen::VectorXd best_sol_ei(dim);

	double bias_weight = psoWeight;
	fill_vector(&best_sol.get(0),best_sol_ei);
	Eigen::VectorXd gb_vec = best_sol_ei-xmean;
	double gb_vec_length = sqrt(gb_vec.transpose() * gb_vec);
	Eigen::VectorXd b_main = B.col(dim-1);
	Eigen::VectorXd bias(dim);
	bias.setZero();

	// Rotation Matrix
	Eigen::MatrixXd R(dim,dim);

	if (gb_vec_length > 0.0)
	{
	    if(sigma < gb_vec_length)
	    {
	    	if(sigma/gb_vec_length <= t_c*gb_vec_length)
	    	{bias = 0.5*gb_vec;}
	    	else
	    	{bias = sigma*gb_vec/gb_vec_length;}
	    }
	    else
	    {bias.setZero();}
	}

	  xmean = xmean + bias;

	  if (psoWeight < 1.0)
	  {
		  Eigen::MatrixXd B_rot(dim,dim);
		  Eigen::DiagonalMatrix<double,Eigen::Dynamic> D_square(dim);

		  Eigen::VectorXd gb_vec_old = best_sol_ei - xold;
		  create_rotmat(b_main,gb_vec_old,R);
		  for (size_t i = 0 ; i < dim ; i++)
		  {B_rot.col(i) = R*B.col(i);}

		  for (size_t i = 0 ; i < dim ; i++)
		  {D_square.diagonal()[i] = D.diagonal()[i] * D.diagonal()[i];}
		  C_pso = B_rot * D_square * B_rot.transpose();

		  Eigen::MatrixXd trUp = C_pso.triangularView<Eigen::Upper>();
		  Eigen::MatrixXd trDw = C_pso.triangularView<Eigen::StrictlyUpper>();
		  C_pso = trUp + trDw.transpose();
	  }
}

void broadcast_best_solution(particle_type & vd,
							 openfpm::vector<double> & best_sol,
							 double & best,
							 double best_sample,
							 openfpm::vector<double> & best_sample_sol)
{
	best_sol.resize(dim);
	auto & v_cl = create_vcluster();

	double best_old = best_sample;
	v_cl.min(best_sample);
	v_cl.execute();

	// The old solution remain the best
	if (best < best_sample)
	{return;}

	best = best_sample;

	size_t rank;
	if (best_old == best_sample)
	{
		rank = v_cl.getProcessUnitID();

		// we own the minimum and we decide who broad cast
		v_cl.min(rank);
		v_cl.execute();

		if (rank == v_cl.getProcessUnitID())
		{
			for (size_t i = 0 ; i < dim ; i++)
			{best_sol.get(i) = best_sample_sol.get(i);}
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

void cmaes_myprctile(openfpm::vector<fun_index> & f_obj, double (& perc)[2], double (& res)[2])
{
	double sar[lambda];
	double availablepercentiles[lambda];
	int idx[hist_size];
	int i,k;

	for (size_t i = 0 ; i < lambda ; i++)
	{
		availablepercentiles[i] = 0.0;
		sar[i] = f_obj.get(i).f;
	}
	std::sort(&sar[0],&sar[lambda]);

	for (size_t i = 0 ; i < 2 ; i++)
	{
		if (perc[i] <= (100.0*0.5/lambda))
		{res[i] = sar[0];}
		else if (perc[i] >= (100.0*(lambda-0.5)/lambda) )
		{res[i] = sar[lambda-1];}
		else
		{
			for (size_t j = 0 ; j < lambda ; j++)
			{availablepercentiles[j] = 100.0 * ((double(j)+1.0)-0.5) / lambda;}

			for (k = 0 ; k < lambda ; k++)
			{if(availablepercentiles[k] >= perc[i]) {break;}}
			k-=1;

			res[i] = sar[k] + (sar[k+1]-sar[k]) * (perc[i]
							-availablepercentiles[k]) / (availablepercentiles[k+1] - availablepercentiles[k]);
		}
	}
}

double maxval(double (& buf)[hist_size], bool (& mask)[hist_size])
{
	double max = 0.0;
	for (size_t i = 0 ; i < hist_size ; i++)
	{
		if (buf[i] > max && mask[i] == true)
		{max = buf[i];}
	}

	return max;
}

double minval(double (& buf)[hist_size], bool (& mask)[hist_size])
{
	double min = std::numeric_limits<double>::max();
	for (size_t i = 0 ; i < hist_size ; i++)
	{
		if (buf[i] < min && mask[i] == true)
		{min = buf[i];}
	}

	return min;
}

void cmaes_intobounds(Eigen::VectorXd & x, Eigen::VectorXd & xout,bool (& idx)[dim], bool & idx_any)
{
	idx_any = false;
	for (size_t i = 0; i < dim ; i++)
	{
		if(x(i) < -5.0)
		{
			xout(i) = -5.0;
			idx[i] = true;
			idx_any = true;
		}
		else if (x(i) > 5.0)
		{
			xout(i) = 5.0;
			idx[i] = true;
			idx_any = true;
		}
		else
		{
			xout(i) = x(i);
			idx[i] = false;
		}
	}
}

void cmaes_handlebounds(openfpm::vector<fun_index> & f_obj,
						double sigma,
						double & validfit,
						Eigen::VectorXd (& arxvalid)[lambda],
						Eigen::VectorXd (& arx)[lambda],
						Eigen::MatrixXd & C,
						Eigen::VectorXd & xmean,
						Eigen::VectorXd & xold,
						double (& weight)[dim],
						double (& fithist)[hist_size],
						bool & iniphase,
						double & validfitval,
						double mu_eff,
						int step,
						int last_restart)
{
	double val[2];
	double value;
	double diag[dim];
	double meandiag;
	int i,k,maxI;
	bool mask[hist_size];
	bool idx[dim];
	Eigen::VectorXd tx(dim);
	int dfitidx[hist_size];
	double dfitsort[hist_size];
	double prct[2] = {25.0,75.0};
	bool idx_any;

	for (size_t i = 0 ; i < hist_size ; i++)
	{
		dfitsort[i] = 0.0;
		dfitidx[i] = 0;

		if (fithist[i] > 0.0)
		{mask[i] = true;}
		else
		{mask[i] = false;}
	}

	for (size_t i = 0 ; i < dim ; i++)
	{diag[i] = C(i,i);}

	maxI = 0;

	meandiag = C.trace()/dim;

	cmaes_myprctile(f_obj, prct, val);
	value = (val[1] - val[0]) / dim / meandiag / (sigma*sigma);

	if (value >= std::numeric_limits<double>::max())
	{
		auto & v_cl = create_vcluster();
		std::cout << "Process " << v_cl.rank() << " warning: Non-finite fitness range" << std::endl;
		value = maxval(fithist,mask);
	}
	else if(value == 0.0)
	{
		value = minval(fithist,mask);
	}
	else if (validfit == 0.0)
	{
		for (size_t i = 0 ; i < hist_size ; i++)
		{fithist[i] = -1.0;}
		validfit = 1;
	}

	for (size_t i = 0; i < hist_size ; i++)
	{
		if(fithist[i] < 0.0)
		{
			fithist[i] = value;
			maxI = i;
			break;
		}
		else if(i == hist_size-1)
		{
			for (size_t k = 0 ; k < hist_size-1 ; k++)
			{fithist[k] = fithist[k+1];}
			fithist[i] = value;
			maxI = i;
		}
	}

	cmaes_intobounds(xmean,tx,idx,idx_any);

	if (iniphase)
	{
		if (idx_any)
		{
			if(maxI == 0)
			{value = fithist[0];}
			else
			{
				openfpm::vector<fun_index> fitsort(maxI+1);
				for (size_t i = 0 ; i <= maxI; i++)
				{
					fitsort.get(i).f = fithist[i];
					fitsort.get(i).id = i;
				}

				fitsort.sort();
				for (size_t k = 0; k <= maxI ; k++)
				{fitsort.get(k).f = fithist[fitsort.get(k).id];}

				if ((maxI+1) % 2 == 0)
				{value = (fitsort.get(maxI/2).f+fitsort.get(maxI/2+1).f)/2.0;}
				else
				{value = fitsort.get(maxI/2).f;}
			}
			for (size_t i = 0 ; i < dim ; i++)
			{
				diag[i] = diag[i]/meandiag;
				weight[i] = 2.0002 * value / diag[i];
			}
			if (validfitval == 1.0 && step-last_restart > 2)
			{
				iniphase = false;
			}
		}
	}

	if(idx_any)
	{
		tx = xmean - tx;
		for(size_t i = 0 ; i < dim ; i++)
		{
			idx[i] = (idx[i] && (fabs(tx(i)) > 3.0*std::max(1.0,sqrt(dim)/mu_eff) * sigma * sqrt(diag[i])));
			idx[i] = (idx[i] && (std::copysign(1.0,tx(i)) == std::copysign(1.0,(xmean(i)-xold(i)))) );
		}
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (idx[i] == true)
			{
				weight[i] = pow(1.2,(std::max(1.0,mu_eff/10.0/dim)))*weight[i];
			}
		}
	}
	double arpenalty[lambda];
	for (size_t i = 0 ; i < lambda ; i++)
	{
		arpenalty[i] = 0.0;
		for (size_t j = 0 ; j < dim ; j++)
		{
			arpenalty[i] += weight[j] * (arxvalid[i](j) - arx[i](j))*(arxvalid[i](j) - arx[i](j));
		}
		f_obj.get(i).f += arpenalty[i];
	}
//	fitness%sel = fitness%raw + bnd%arpenalty;
}

double adjust_sigma(double sigma, Eigen::MatrixXd & C)
{
	for (size_t i = 0 ; i < dim ; i++)
	{
		if (sigma*sqrt(C(i,i)) > 5.0)
		{sigma = 5.0/sqrt(C(i,i));}
	}

	return sigma;
}


void cma_step(particle_type & vd, int step,  double & best,
			  int & best_i, openfpm::vector<double> & best_sol,
			  size_t & fun_eval)
{
	size_t fe = 0;
	Eigen::VectorXd xmean(dim);
	Eigen::VectorXd mean_z(dim);
	Eigen::VectorXd arxvalid[lambda];
	Eigen::VectorXd arx[lambda];

	for (size_t i = 0 ; i < lambda ; i++)
	{
		arx[i].resize(dim);
		arxvalid[i].resize(dim);
	}

	double best_sample = std::numeric_limits<double>::max();
	openfpm::vector<double> best_sample_sol(dim);

	openfpm::vector<fun_index> f_obj(lambda);

	int counteval = step*lambda;

	auto it = vd.getDomainIterator();
	while (it.isNext())
	{
		auto p = it.get();

		if (vd.getProp<stop>(p) == true)
		{++it;continue;}

		Eigen::VectorXd (& arz)[lambda] = vd.getProp<Zeta>(p);

		// fill the mean vector;

		fill_vector(vd.getPos(p),xmean);

		for (size_t j = 0 ; j < lambda ; j++)
		{
			vd.getProp<Zeta>(p)[j] = generateGaussianVector<dim>();
			arx[j] = xmean + vd.getProp<sigma>(p)*vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<Zeta>(p)[j];

			// sample point has to be inside -5.0 and 5.0
			for (size_t i = 0 ; i < dim ; i++)
			{
				if (arx[j](i) < -5.0)
				{arxvalid[j](i) = -5.0;}
				else if (arx[j](i) > 5.0)
				{arxvalid[j](i) = 5.0;}
				else
				{arxvalid[j](i) = arx[j](i);}
			}

			f_obj.get(j).f = hybrid_composition<dim>(arxvalid[j]);
			f_obj.get(j).id = j;
			fe++;

			// Get the best ever
			if (f_obj.get(j).f < best_sample)
			{
				best_sample = f_obj.get(j).f;

			    // Copy the new mean as position of the particle
			    for (size_t i = 0 ; i < dim ; i++)
			    {best_sample_sol.get(i) = arxvalid[j](i);}
			}
		}

		// Add penalities for out of bound points
		cmaes_handlebounds(f_obj,vd.getProp<sigma>(p),
						   vd.getProp<validfit>(p),arxvalid,
						   arx,vd.getProp<Cov_m>(p),
						   xmean,vd.getProp<xold>(p),vd.getProp<weight>(p),
						   vd.getProp<fithist>(p),vd.getProp<iniphase>(p),
						   vd.getProp<validfit>(p),mu_eff,
						   step,vd.getProp<last_restart>(p));

		f_obj.sort();

		for (size_t j = 0 ; j < lambda ; j++)
		{vd.getProp<ord>(p)[j] = f_obj.get(j).id;}

		vd.getProp<xold>(p) = xmean;

		// Calculate weighted mean

		xmean.setZero();
		mean_z.setZero();
		for (size_t j = 0 ; j < mu ; j++)
		{
			xmean += weight_sample(j)*arx[vd.getProp<ord>(p)[j]];
			mean_z += weight_sample(j)*vd.getProp<Zeta>(p)[vd.getProp<ord>(p)[j]];
		}

		vd.getProp<xmean_st>(p) = xmean;
		vd.getProp<meanz_st>(p) = mean_z;

		++it;
	}

	// Find the best point across processors
	broadcast_best_solution(vd,best_sol,best,best_sample,best_sample_sol);

	// bool calculate B and D
	bool calc_bd = counteval - eigeneval > lambda/(ccov)/dim/10;
	if (calc_bd == true)
	{eigeneval = counteval;}

	auto it2 = vd.getDomainIterator();
	while (it2.isNext())
	{
		auto p = it2.get();

		if (vd.getProp<stop>(p) == true)
		{++it2;continue;}

		xmean = vd.getProp<xmean_st>(p);
		mean_z = vd.getProp<meanz_st>(p);

		vd.getProp<path_s>(p) = vd.getProp<path_s>(p)*(1.0 - cs) + sqrt(cs*(2.0-cs)*mu_eff)*vd.getProp<B>(p)*mean_z;

		double hsig = vd.getProp<path_s>(p).norm()/sqrt(1.0-pow((1.0-cs),(2.0*double((step-vd.getProp<last_restart>(p))))))/chiN < 1.4 + 2.0/(dim+1);

		vd.getProp<path_c>(p) = (1-cc)*vd.getProp<path_c>(p) + hsig * sqrt(cc*(2-cc)*mu_eff)*(vd.getProp<B>(p)*vd.getProp<D>(p)*mean_z);

		if (step % N_pso == 0)
		{
			Eigen::MatrixXd C_pso(dim,dim);
			updatePso(best_sol,vd.getProp<sigma>(p),xmean,vd.getProp<xold>(p),vd.getProp<B>(p),vd.getProp<D>(p),C_pso);

			// Adapt covariance matrix C
			vd.getProp<Cov_m>(p) = (1.0-ccov+(1.0-hsig)*ccov*cc*(2.0-cc)/mu_eff)*vd.getProp<Cov_m>(p) +
									ccov*(1.0/mu_eff)*(vd.getProp<path_c>(p)*vd.getProp<path_c>(p).transpose());

			for (size_t i = 0 ; i < mu ; i++)
			{vd.getProp<Cov_m>(p) += ccov*(1.0-1.0/mu_eff)*(vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<Zeta>(p)[vd.getProp<ord>(p)[i]])*weight_sample(i)*
										  (vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<Zeta>(p)[vd.getProp<ord>(p)[i]]).transpose();
			}

	    	vd.getProp<Cov_m>(p) = psoWeight*vd.getProp<Cov_m>(p) + (1.0 - psoWeight)*C_pso;
	    }
	    else
	    {
			// Adapt covariance matrix C
			vd.getProp<Cov_m>(p) = (1.0-ccov+(1.0-hsig)*ccov*cc*(2.0-cc)/mu_eff)*vd.getProp<Cov_m>(p) +
									ccov*(1.0/mu_eff)*(vd.getProp<path_c>(p)*vd.getProp<path_c>(p).transpose());

			for (size_t i = 0 ; i < mu ; i++)
			{vd.getProp<Cov_m>(p) += ccov*(1.0-1.0/mu_eff)*(vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<Zeta>(p)[vd.getProp<ord>(p)[i]])*weight_sample(i)*
				                          (vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<Zeta>(p)[vd.getProp<ord>(p)[i]]).transpose();
			}
	    }

		// Numeric error

		double smaller = std::numeric_limits<double>::max();
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (vd.getProp<sigma>(p)*sqrt(vd.getProp<D>(p).diagonal()[i]) > 5.0)
			{
				if (smaller > 5.0/sqrt(vd.getProp<D>(p).diagonal()[i]))
				{smaller = 5.0/sqrt(vd.getProp<D>(p).diagonal()[i]);}
			}
		}
		if (smaller != std::numeric_limits<double>::max())
		{vd.getProp<sigma>(p) = smaller;}

		//Adapt step-size sigma
		vd.getProp<sigma>(p) = vd.getProp<sigma>(p)*exp((cs/d_amps)*(vd.getProp<path_s>(p).norm()/chiN - 1));

		// Update B and D from C

		if (calc_bd)
		{
			Eigen::MatrixXd trUp = vd.getProp<Cov_m>(p).triangularView<Eigen::Upper>();
			Eigen::MatrixXd trDw = vd.getProp<Cov_m>(p).triangularView<Eigen::StrictlyUpper>();
			vd.getProp<Cov_m>(p) = trUp + trDw.transpose();

			// Eigen decomposition
			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_solver;

			eig_solver.compute(vd.getProp<Cov_m>(p));

			for (size_t i = 0 ; i < eig_solver.eigenvalues().size() ; i++)
			{vd.getProp<D>(p).diagonal()[i] = sqrt(eig_solver.eigenvalues()[i]);}
			vd.getProp<B>(p) = eig_solver.eigenvectors();

			// Make first component always positive
			for (size_t i = 0 ; i < dim ; i++)
			{
				if (vd.getProp<B>(p)(0,i) < 0)
				{vd.getProp<B>(p).col(i) = - vd.getProp<B>(p).col(i);}
			}

			Eigen::MatrixXd tmp = vd.getProp<B>(p).transpose();
		}

	    // Copy the new mean as position of the particle
	    for (size_t i = 0 ; i < dim ; i++)
	    {vd.getPos(p)[i] = xmean(i);}

	    vd.getProp<sigma>(p) = adjust_sigma(vd.getProp<sigma>(p),vd.getProp<Cov_m>(p));

	    // Stop conditions
	    bool stop_tol = true;
	    bool stop_tolX = true;
	    for (size_t i = 0 ; i < dim ; i++)
	    {
	    	stop_tol &= (vd.getProp<sigma>(p)*std::max(fabs(vd.getProp<path_c>(p)(i)),sqrt(vd.getProp<Cov_m>(p)(i,i)))) < stopTolX;
	    	stop_tolX &= vd.getProp<sigma>(p)*sqrt(vd.getProp<D>(p).diagonal()[i]) > stopTolUpX;
	    }

	    vd.getProp<stop>(p) = stop_tol | stop_tolX;

	    // Escape flat fitness, or better terminate?
	    if (f_obj.get(0).f == f_obj.get(std::ceil(0.7*lambda)).f )
	    {
	    	vd.getProp<sigma>(p) = vd.getProp<sigma>(p)*exp(0.2+cs/d_amps);
	    	std::cout << "warning: flat fitness, consider reformulating the objective";

	    	// Stop it
	    	vd.getProp<stop>(p) = true;
	    }

	    if (vd.getProp<stop>(p) == true)
	    {std::cout << "Stopped" << std::endl;}

	    if (restart_cma && vd.getProp<stop>(p) == true)
	    {
	    	std::cout << "------- Restart #" << std::endl;

	    	std::cout << "---------------------------------" << std::endl;
	    	std::cout << "Best: " << best << "   " << fun_eval << std::endl;
	        std::cout << "---------------------------------" << std::endl;

	        vd.getProp<last_restart>(p) = step;
	        vd.getProp<xold>(p).setZero();

			for (size_t i = 0 ; i < vd.getProp<D>(p).diagonal().size() ; i++)
			{vd.getProp<D>(p).diagonal()[i] = 1.0;}
			vd.getProp<B>(p).resize(dim,dim);
			vd.getProp<B>(p).setIdentity();
			vd.getProp<Cov_m>(p) = vd.getProp<B>(p)*vd.getProp<D>(p)*vd.getProp<D>(p)*vd.getProp<B>(p);
			vd.getProp<path_s>(p).resize(dim);
			vd.getProp<path_s>(p).setZero(dim);
			vd.getProp<path_c>(p).resize(dim);
			vd.getProp<path_c>(p).setZero(dim);
			vd.getProp<stop>(p) = false;
			vd.getProp<iniphase>(p) = true;
			vd.getProp<last_restart>(p) = 0;
			vd.getProp<sigma>(p) = 2.0;

			// a different point in space
			for (size_t i = 0 ; i < dim ; i++)
			{
				// we define x, assign a random position between 0.0 and 1.0
				vd.getPos(p)[i] = 10.0*(double)rand() / RAND_MAX - 5.0;
			}

			// Initialize the bound history

			for (size_t i = 0 ; i < hist_size ; i++)
			{vd.getProp<fithist>(p)[i] = -1.0;}
			vd.getProp<fithist>(p)[0] = 1.0;
			vd.getProp<validfit>(p) = 0.0;
		}

		++it2;
	}

	auto & v_cl = create_vcluster();
	v_cl.sum(fe);
	v_cl.execute();

	fun_eval += fe;
}



int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

	auto & v_cl = create_vcluster();

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

	prepare_f15<dim>();

	// extended boundary around the domain, and the processor domain
	Ghost<dim,double> g(0.0);

    particle_type vd(16,domain,bc,g);

    // Initialize constants

    stop_fitness = 1e-10;
    size_t stopeval = 1e3*dim*dim;

    // Strategy parameter setting: Selection
    init_weight();

    // Strategy parameter setting: Adaptation
    cc = 4.0 / (dim+4.0);
    cs = (mu_eff+2.0) / (double(dim)+mu_eff+3.0);
    ccov = (1.0/mu_eff) * 2.0/((dim+1.41)*(dim+1.41)) +
    	   (1.0 - 1.0/mu_eff)* std::min(1.0,(2.0*mu_eff-1.0)/((dim+2.0)*(dim+2.0) + mu_eff));
    d_amps = 1 + 2*std::max(0.0, sqrt((mu_eff-1.0)/(dim+1))-1) + cs;

    chiN = sqrt(dim)*(1.0-1.0/(4.0*dim)+1.0/(21.0*dim*dim));

	//! \cond [assign position] \endcond


	// initialize the srand
	int seed = 24756*v_cl.rank()*v_cl.rank() + time(NULL);
	srand(seed);

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		for (size_t i = 0 ; i < dim ; i++)
		{
			// we define x, assign a random position between 0.0 and 1.0
			vd.getPos(p)[i] = 10.0*(double)rand() / RAND_MAX - 5.0;
		}

		vd.getProp<sigma>(p) = 2.0;

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
		vd.getProp<stop>(p) = false;
		vd.getProp<iniphase>(p) = true;
		vd.getProp<last_restart>(p) = 0;

		// Initialize the bound history

		for (size_t i = 0 ; i < hist_size ; i++)
		{vd.getProp<fithist>(p)[i] = -1.0;}
		vd.getProp<fithist>(p)[0] = 1.0;
		vd.getProp<validfit>(p) = 0.0;

		// next particle
		++it;
	}

	if (v_cl.rank() == 0)
	{std::cout << "Starting PS-CMA-ES" << std::endl;}

	double best = 0.0;
	int best_i = 0;

	best = std::numeric_limits<double>::max();
	openfpm::vector<double> best_sol(dim);
	// now do several iteration

	int stop_cond = 0;
	size_t fun_eval = 0;
	int i = 0;
	while (fun_eval < max_fun_eval && best > 120.000001)
	{
		// sample offspring
		cma_step(vd,i+1,best,best_i,best_sol,fun_eval);

		i++;
	}

	if (v_cl.rank() == 0)
	{
		std::cout << "Best solution: " << best << " with " << fun_eval << std::endl;
		std::cout << "at: " << std::endl;

		for (size_t i = 0 ; i < best_sol.size() ; i++)
		{
			std::cout << best_sol.get(i) << " ";
		}
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
