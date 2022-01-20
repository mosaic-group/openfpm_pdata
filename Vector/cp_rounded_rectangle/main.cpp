#include <math.h>
#include <sys/_types/_size_t.h>
#include "Vector/vector_dist.hpp"
#include "DCPSE/Dcpse.hpp"
#include "../../openfpm_numerics/src/DCPSE/MonomialBasis.hpp"
#include "Draw/DrawParticles.hpp"
#include "../../openfpm_numerics/src/level_set/particle_cp/particle_cp.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

const double dp = 1/2048.0;
const double A = 0.5;	  //length of the rectangle part
const double B = 0.25;	  //radius of the rounded part
const double omega = 2.0; //parameter for steering the errors in the initialization based on sines for the rectangle part
const double band_width = 10.0*dp;

const int sdf = 0;
const int sdfgrad = 1;
const int curvature = 2;
const int surf_flag = 3;
const int num_neibs = 4;
const int min_sdf = 5;
const int min_sdf_x = 6;
const int interpol_coeff = 7;
const int sdf_analytical = 8;

typedef vector_dist<2, double, aggregate<double, double[2], double, int, int, double, double[2], double[16], double>> particles;
//										|		|	 	 |         |    |		 |		|	|						|
//									     sdf  sdfgrad curvature surf_flag num_neibs min sdf in support	min sdf location poly coefficients analytical sdf

//typedef vector_dist<2, double, aggregate<vect_dist_key_dx>> particles_surface;

inline double get_distance(double x, double y)
{
	double dist;

	x = abs(x);
	y = abs(y);
	
	if(x<A/2) // point is in the rectangle part
	{
		dist = B - y;
	}
	else
	{
		dist = B - sqrt((x - A/2)*(x - A/2) + y*y);
	}
	return(dist);
}

inline double guess_sdf(double x, double y)
{
	double sdf_guess;
	x = abs(x);
	y = abs(y);

	if(x<A/2)
	{
		//sdf_guess = - sin(omega*(y - A/2));
		sdf_guess = -(y*y - B*B);
	}
	else
	{
		sdf_guess = -((x - A/2)*(x - A/2) + y*y - B*B);
	}

	return(sdf_guess);
}

int return_sign(double phi)
{
	if (phi > 0) return 1;
	if (phi < 0) return -1;
	return 0;
}

double randZeroToOne()
{
    return (rand() / (RAND_MAX + 1.));
}

inline void perturb_pos(particles & vd)
{
	auto part = vd.getDomainIterator();

	while(part.isNext())
	{
		auto a = part.get();

		const double x = vd.getPos(a)[0];
		const double y = vd.getPos(a)[1];

		vd.getPos(a)[0] += randZeroToOne()*dp/3.0;
		vd.getPos(a)[1] += randZeroToOne()*dp/3.0;

		++part;
	}
}

template <typename CellList> inline void update_sdfs(particles & vd, CellList & NN)//
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
	{
		auto a = part.get();
		Point<2,double> xa = vd.getPos(a);
		double x = xa[0];
		double y = xa[1];
		
		vd.getProp<sdf_analytical>(a) = get_distance(x,y);
		vd.getProp<sdf>(a) = guess_sdf(x,y);
		++part;
	}
}

//################################################################# main function ##################
int main(int argc, char* argv[])
{
	openfpm_init(&argc, &argv);

	const double l = 4.0;
	Box<2, double> domain({-l/2.0, -l/2.0}, {l/2.0, l/2.0});
	size_t sz[2] = {(size_t)(l/dp + 0.5), (size_t)(l/dp + 0.5)};

	size_t bc[2] = {NON_PERIODIC, NON_PERIODIC};
	Ghost<2, double> g(0.0);

	particles vd(0, domain, bc, g, DEC_GRAN(512));
	//particles_surface vd_s(vd.getDecomposition(), 0);

	openfpm::vector<std::string> names({"sdf", "sdfgradient", "curvature", "surface_flag", "number_neibs", "min_sdf_neibors", "min_sdf_location", "interpol_coeff", "sdf_analytical"});
	vd.setPropNames(names);

	Box<2, double> particle_box({-l/2.0, -l/2.0}, {l/2.0, l/2.0});
	auto particle_it = DrawParticles::DrawBox(vd, sz, domain, particle_box);

	while (particle_it.isNext())
	{
		double x = particle_it.get().get(0);
		double y = particle_it.get().get(1);

		double dist = get_distance(x, y);

		//std::cout<<x<<"\t"<<y<<"\t"<<dist<<std::endl;

		if (abs(dist) < band_width/2.0)
		{
			vd.add();
			vd.getLastPos()[0] = x;
			vd.getLastPos()[1] = y;
			double sdf_guess = guess_sdf(x, y);
			vd.template getLastProp<sdf>() = sdf_guess;
			vd.template getLastProp<sdf_analytical>() = dist;
			vd.template getLastProp<surf_flag>() = 0;
			vd.template getLastProp<num_neibs>() = 1;

		}

		++particle_it;

	}

	auto NN = vd.getCellList(2.6*dp);
	
	perturb_pos(vd);
	vd.map();
	NN = vd.getCellList(2.6*dp);
	update_sdfs(vd, NN);
	vd.write("rounded_rectangle_init");

	Redist_options rdistoptions;
	rdistoptions.max_iter = 1000;
	rdistoptions.tolerance = 1e-13;
	rdistoptions.H = dp;
	rdistoptions.r_cutoff_factor = 2.5;
	rdistoptions.sampling_radius = 0.75*band_width;
	rdistoptions.polynomial_degree = "bicubic";
	//rdistoptions.barrier_coefficient = -5*1e-6;//-8.75*1e-7;//-10
	//rdistoptions.armijo_tau = 0.9;
	//rdistoptions.armijo_c = 0.9;
	//rdistoptions.support_prevent = 0.9;
	//rdistoptions.fail_projection = 1;
	//rdistoptions.init_project = 1;
	particle_cp_redistancing pcprdist(vd, rdistoptions);

	pcprdist.run_redistancing();
	vd.write("rounded_rectangle_after_redistancing");

	openfpm_finalize();
}
