#define SE_CLASS1

#include <math.h>
#include <sys/_types/_size_t.h>
#include "Vector/vector_dist.hpp"
#include "DCPSE/Dcpse.hpp"
#include "DCPSE/MonomialBasis.hpp"
#include "Draw/DrawParticles.hpp"
#include "../../openfpm_numerics/src/level_set/particle_cp/particle_cp.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"
#include <chrono>

const int sdf = 0;
const int normal = 1;
const int curvature = 2;
const int surf_flag = 3;
const int rho = 4;
const int p = 5;
const int force = 6;
const int vel = 7;
const int vel_prev = 8;

typedef vector_dist<3, double, aggregate<double, double[3], double, int, double, double, double[3], double[3], double[3]>> particles;
//											|		|	 	 |         |    |		 |		|	|						|
//										     sdf  sdfgrad curvature surf_flag 	density, previous density	pressure viscosity (dyn), force, velocity, previous velocity

//typedef vector_dist<2, double, aggregate<vect_dist_key_dx>> particles_surface;

const double dp = 1/32.0;
const double H = std::sqrt(3.0*dp*dp);
const double gamma_eos = 7.0;
const double c = 20.0;
const double rho_0 = 1.0;
const double rho_1 = 1.0;
const double eta_phase0 = 0.2;
const double eta_phase1 = 0.2;
const double alpha = 1.0;
const double p_background = 0.0;

const double band_width = 6.0*dp;

const double l_0 = 0.6;
const double l_1 = 1.0;
const double t_end = 1.0;
const double M_0 = std::pow(l_0, 3)*rho_0;
const double M_1 = (std::pow(l_1, 3) - std::pow(l_0, 3))*rho_1; 
const double np_0 = std::pow(l_0/dp, 3);
const double np_1 = std::pow(l_1/dp, 3) - np_0;
const double m_0 = M_0/np_0;
const double m_1 = M_1/np_1;

int return_sign(double phi)
{
	if (phi > 0) return 1;
	if (phi < 0) return -1;
	return 0;
}

inline double getDistancePointCube(double x, double y, double z)
{
	double d;
	if (x<=l_0/2.0)
	{
	     if (y<=l_0/2.0)
	          d=z-l_0/2.0;
	     else
	     {
	          if (z<=l_0/2.0)
	                d=y-l_0/2.0;
	          else
	                d=std::sqrt(std::pow(y-l_0/2.0, 2.0)+std::pow(z-l_0/2.0, 2.0));
	     }
	}
	else
	{
	     if (y<=l_0/2.0)
	     {
	          if (z<=l_0/2.0)
	               d=x - l_0/2.0;
	          else
	               d=std::sqrt(std::pow(x-l_0/2.0, 2.0)+std::pow(z-l_0/2.0, 2.0));
	     }
	     else
	          if (z<=l_0/2.0)
	               d=std::sqrt(std::pow(x-l_0/2.0, 2.0)+std::pow(y-l_0/2.0, 2.0));
	          else
	               d=std::sqrt(std::pow(x-l_0/2.0, 2.0)+std::pow(y-l_0/2.0, 2.0)+std::pow(z-l_0/2.0, 2.0));
	     }
	return d;
}

// kernel function

const double a2 = 1.0/M_PI/H/H/H;

inline double Wab(double r)
{
	r /= H;

	if (r < 1.0)
		return (1.0 - 3.0/2.0*r*r + 3.0/4.0*r*r*r)*a2;
	else if (r < 2.0)
		return (1.0/4.0*(2.0 - r*r)*(2.0 - r*r)*(2.0 - r*r))*a2;
	else
		return 0.0;
}

const double c1 = -3.0/M_PI/H/H/H/H;
const double d1 = 9.0/4.0/M_PI/H/H/H/H;
const double c2 = -3.0/4.0/M_PI/H/H/H/H;
const double a2_4 = 0.25*a2;
// Filled later
double W_dap = 0.0;

inline void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print, double & dwdrab)
{
	const double qq=r/H;

    double qq2 = qq * qq;
    double fac1 = (c1*qq + d1*qq2)/r;
    double b1 = (qq < 1.0)?1.0f:0.0f;

    double wqq = (2.0 - qq);
    double fac2 = c2 * wqq * wqq / r;
    double b2 = (qq >= 1.0 && qq < 2.0)?1.0f:0.0f;

    double factor = (b1*fac1 + b2*fac2);

    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);

    dwdrab = factor;
}

template <typename CellList> inline void density_summation(particles & vd, CellList & NN)
{
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
	{
		auto a = part.get();
		Point<3, double> xa = vd.getPos(a);
		double rho_a = 0.0;
		
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
		
		while (Np.isNext() == true)
		{
			auto b = Np.get();
			Point<3, double> xb = vd.getPos(b);
			Point<3, double> dr = xa - xb;
			double r2 = norm2(dr);
			
			if (r2 < 4.0*H*H)
			{
				double r = std::sqrt(r2);
				rho_a += Wab(r); 
			}
			++Np;
		}
		vd.getProp<rho>(a) = m_0*rho_a;
		++part;
	}
}

inline double EqOfState(double dens, int phase)
{
	double ref_dens = 0.0;
	if (phase < 0) ref_dens = rho_0;
	else if (phase > 0) ref_dens = rho_1;
	double p;
	p = (c*c*ref_dens/gamma_eos)*(std::pow(dens/ref_dens, gamma_eos) - 1) + p_background;
	
	return(p);
}

template<typename CellList> inline void calc_forces(particles & vd, CellList & NN, double & max_visc)
{
	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// Get the position xp of the particle
		Point<3,double> xa = vd.getPos(a);

		// Take the mass of the particle dependently if it is FLUID 0 or 1
		double massa = (return_sign(vd.getProp<sdf>(a)) < 0)?m_0:m_1;
		double etaa = (return_sign(vd.getProp<sdf>(a)) < 0)?eta_phase0:eta_phase1;
		// Get the density of the of the particle a
		double rhoa = vd.getProp<rho>(a);
		double Va = massa/rhoa;

		// Get the pressure of the particle a
		double Pa = vd.getProp<p>(a);

		// Get the Velocity of the particle a
		Point<3,double> va = vd.getProp<vel>(a);

		// Get an iterator over the neighborhood particles of p
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

		// For each neighborhood particle
		while (Np.isNext() == true)
		{
			// ... q
			auto b = Np.get();

			// Get the position xp of the particle
			Point<3,double> xb = vd.getPos(b);

			// if (p == q) skip this particle
			if (a.getKey() == b)	{++Np; continue;};

			double massb = (return_sign(vd.getProp<sdf>(a)) < 0)?m_0:m_1;
			double etab = (return_sign(vd.getProp<sdf>(a)) < 0)?eta_phase0:eta_phase1;
			Point<3,double> vb = vd.getProp<vel>(b);
			double Pb = vd.getProp<p>(b);
			double rhob = vd.getProp<rho>(b);
			double Vb = massb/rhob;

			// Get the distance between p and q
			Point<3,double> dr = xa - xb;
			// take the norm of this vector
			double r2 = norm2(dr);

			// if they interact
			if (r2 < 4.0*H*H)
			{
				double r = sqrt(r2);

				Point<3,double> v_rel = va - vb;

				Point<3,double> DW;
				double dwdrab;
				DWab(dr,DW,r,false,dwdrab);

				double factor_p = - (Va*Va + Vb*Vb)*(rhoa*Pb + rhob*Pa)/(rhoa + rhob);
				double factor_visc = (2*etaa*etab/(etaa + etab))*(Va*Va + Vb*Vb)*dwdrab/r;

				vd.getProp<force>(a)[0] += factor_p * DW.get(0) + factor_visc*v_rel[0];
				vd.getProp<force>(a)[1] += factor_p * DW.get(1) + factor_visc*v_rel[1];
				vd.getProp<force>(a)[2] += factor_p * DW.get(2) + factor_visc*v_rel[2];
			}

			++Np;
		}

		vd.getProp<force>(a)[0] = vd.getProp<force>(a)[0]/massa;
		vd.getProp<force>(a)[1] = vd.getProp<force>(a)[1]/massa;
		vd.getProp<force>(a)[2] = vd.getProp<force>(a)[2]/massa;

		if (std::abs(vd.getProp<sdf>(a)) < (2.0*H))
		{
			double smoothing_factor = Wab(std::abs(vd.getProp<sdf>(a)))/rhoa; //latter is part of the equations, not the smoothing factor itself, but for convinience its in there
			vd.getProp<force>(a)[0] = vd.getProp<force>(a)[0] - smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[0];
			vd.getProp<force>(a)[1] = vd.getProp<force>(a)[1] - smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[1];
			vd.getProp<force>(a)[2] = vd.getProp<force>(a)[2] - smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[2];
		}

		++part;
	}
}

int main(int argc, char* argv[])
{
	openfpm_init(&argc, &argv);

	Box<3, double> domain({-l_1/2.0, -l_1/2.0, -l_1/2.0}, {l_1/2.0, l_1/2.0, l_1/2.0});
	size_t sz[3] = {(size_t)(l_1/dp + 0.5),(size_t)(l_1/dp + 0.5),(size_t)(l_1/dp + 0.5)};

	size_t bc[3] = {PERIODIC, PERIODIC, PERIODIC};
	Ghost<3, double> g(band_width);

	particles vd(0, domain, bc, g, DEC_GRAN(512));
	//particles_surface vd_s(vd.getDecomposition(), 0);

	openfpm::vector<std::string> names({"sdf","normal","curvature","surf_flag","rho","rho_prev","p","eta","f","vel","vel_prev"});
	vd.setPropNames(names);
	Box<3, double> particle_box({-l_1/2.0 + dp/2.0, -l_1/2.0 + dp/2.0, -l_1/2.0 + dp/2.0}, {l_1/2.0 - dp/2.0, l_1/2.0 - dp/2.0, l_1/2.0 - dp/2.0});


	auto particle_it = DrawParticles::DrawBox(vd, sz, domain, particle_box);

	while (particle_it.isNext())
	{
		double x = particle_it.get().get(0);
		double y = particle_it.get().get(1);
		double z = particle_it.get().get(2);

		vd.add();
		vd.getLastPos()[0] = x;
		vd.getLastPos()[1] = y;
		vd.getLastPos()[2] = z;
		vd.template getLastProp<surf_flag>() = 0;

		double dist = getDistancePointCube(abs(x), abs(y), abs(z));
		//if ((abs(x) < l_0/2.0) && (abs(y) < l_0/2.0) && (abs(z) < l_0/2.0)) dist = -1.0*dist; // we're inside then

		vd.template getLastProp<sdf>() = dist;
		vd.template getLastProp<normal>()[0] = 0.0;
		vd.template getLastProp<normal>()[1] = 0.0;
		vd.template getLastProp<normal>()[2] = 0.0;
		vd.template getLastProp<curvature>() = 0.0;
		vd.template getLastProp<vel>()[0] = 0.0;
		vd.template getLastProp<vel>()[1] = 0.0;
		vd.template getLastProp<vel>()[2] = 0.0;
		vd.template getLastProp<vel_prev>()[0] = 0.0;
		vd.template getLastProp<vel_prev>()[1] = 0.0;
		vd.template getLastProp<vel_prev>()[2] = 0.0;
		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;
		vd.template getLastProp<rho>() = 0.0;
		vd.template getLastProp<p>() = 0.0;

		++particle_it;
	}

	vd.map();
	vd.write("droplet_init");
}




