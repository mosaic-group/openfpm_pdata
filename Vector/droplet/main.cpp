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
const int pos_prev = 9;

typedef vector_dist<3, double, aggregate<double, double[3], double, int, double, double, Point<3, double>, Point<3, double>, Point<3, double>, Point<3, double>>> particles;
//											|		|	 	 |         |    |		 |		|	|						|
//										     sdf  sdfgrad curvature surf_flag 	density, previous density	pressure viscosity (dyn), force, velocity, previous velocity

//typedef vector_dist<2, double, aggregate<vect_dist_key_dx>> particles_surface;

const double dp = 1/32.0;
//const double H = std::sqrt(3.0*dp*dp);
const double H = 1.3*dp;
const double gamma_eos = 1.4;
const double c = 200.0;
const double rho_0 = 1.0;
const double rho_1 = 1.0;
const double eta_phase0 = 0.2;
const double eta_phase1 = 0.2;
const double alpha = 10.0;
const double p_background = 0.0;

const double dt = 1.0*1e-4;
double t;
int corr;
const int colorfield = 0;

const double band_width = 10.0*dp;

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
	if ((x<=l_0/2.0) && (y<=l_0/2.0) && (z<=l_0/2.0))
	{
		d = std::min(std::min(std::abs(x-l_0/2.0), std::abs(y-l_0/2.0)), std::abs(z-l_0/2.0));
	}
	return d;
}

// kernel function

const double a2 = 21.0/(16.0*M_PI*H*H*H);
inline double Wab(double r)
{
	const double q = r/H;
	if (q <= 2.0)
	{
		double factor = 1.0 - q/2.0;
		factor = factor*factor;
		factor = factor*factor;
		return(a2*factor*(1.0 + 2.0*q));
	}
	else return(0.0);
}

const double a3 = 5.0/(8.0*H);
inline double Wab1D(double r)
{
    const double q = r/H;
    if (q <= 2.0)
    {
        double factor = 1.0 - q/2.0;
        factor = factor*factor*factor;
        return(a3*factor*(1.5*q + 1));
    }
    else return(0.0);
}

inline void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print, double & dwdrab)
{
	const double qq=r/H;
	if (qq <= 2.0)
    {
		double factor = (-5.0*a2/(H))*qq*(1.0 - qq/2.0)*(1.0 - qq/2.0)*(1.0 - qq/2.0);

    	DW.get(0) = factor * dx.get(0)/r;
    	DW.get(1) = factor * dx.get(1)/r;
    	DW.get(2) = factor * dx.get(2)/r;

    	dwdrab = factor;
    }
	else
	{
		DW.get(0) = 0.0;
		DW.get(1) = 0.0;
		DW.get(2) = 0.0;

		dwdrab = 0.0;
	}
}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPH geometric computing part %%%%%%%%%%%%%%%%%%%%%
const double eps_normal = 0.01/H;
inline double coloraverage(double rhoa, double rhob, int typea, int typeb)
{
    double cij;
    if (typea==typeb)
    {cij = 0.0;}
    else cij = 1.0;
    return (rhoa/(rhoa+rhob)*cij);
}

template<typename CellList> inline void calc_surface_normals(particles & vd, CellList & NN) {
    vd.template ghost_get<rho, sdf>();
    auto part = vd.getDomainIterator();

    // Update the cell-list
    vd.updateCellList(NN);

    // For each particle ...
    while (part.isNext()) {
        // ... a
        auto a = part.get();
        vd.getProp<normal>(a)[0] = 0.0;
        vd.getProp<normal>(a)[1] = 0.0;
        //if (std::abs(vd.getProp<sdf>(a))>2.0*H) {++part; continue;};
        // Get the position xp of the particle
        Point<2, double> xa = vd.getPos(a);

        double massa = (return_sign(vd.getProp<sdf>(a)) < 0) ? m_0 : m_1;
        double ca = return_sign(vd.getProp<sdf>(a));
        // Get the density of the particle a
        double rhoa = vd.getProp<rho>(a);

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

        // For each neighborhood particle
        while (Np.isNext() == true) {
            // ... q
            auto b = Np.get();

            // Get the position xp of the particle
            Point<2, double> xb = vd.getPos(b);

            // if (p == q) skip this particle

            if (a.getKey() == b) {
                ++Np;
                continue;
            };

            double massb = (return_sign(vd.getProp<sdf>(b)) < 0) ? m_0 : m_1;

            double rhob = vd.getProp<rho>(b);

            // Get the distance between p and q
            Point<2, double> dr = xa - xb;
            // take the norm of this vector
            double r2 = norm2(dr);

            // if they interact
            if (r2 < 4.0 * H * H) {
                double r = sqrt(r2); // norm2 is norm^2

                Point<2, double> DW;
                double dwdrab;
                DWab(dr, DW, r, false, dwdrab);
                double cfactor = 0.0;
                double cb = return_sign(vd.getProp<sdf>(b));
                if (colorfield) {
                    //cfactor = rhoa/massa*((massa/rhoa)*(massa/rhoa)+(massb/rhob)*(massb/rhob))*coloraverage(rhoa,rhob,return_sign(vd.getProp<sdf>(a)),return_sign(vd.getProp<sdf>(b)));
                    cfactor = massb / rhob * (cb - ca);
                } else cfactor = massb / rhob * (vd.getProp<sdf>(b) - vd.getProp<sdf>(a));
                vd.getProp<normal>(a)[0] += cfactor * DW.get(0);
                vd.getProp<normal>(a)[1] += cfactor * DW.get(1);
            }
            ++Np;
        }
        //std::cout<<"normal:\n"<<vd.getProp<normal>(a)[0]<<", "<<vd.getProp<normal>(a)[1]<<std::endl;
        // normalize normal to obtain unit surface normal
        if (false) {
            double colornorm = sqrt(vd.getProp<normal>(a)[0] * vd.getProp<normal>(a)[0] +
                                    vd.getProp<normal>(a)[1] * vd.getProp<normal>(a)[1]);
            if (colornorm < eps_normal) {
                vd.getProp<normal>(a)[0] = 0.0;
                vd.getProp<normal>(a)[0] = 0.0;
            } else {
                vd.getProp<normal>(a)[0] = vd.getProp<normal>(a)[0] / colornorm;
                vd.getProp<normal>(a)[1] = vd.getProp<normal>(a)[1] / colornorm;
            }
        }

        if (colorfield) {
            if (return_sign(vd.getProp<sdf>(a)) > 0) {
                //vd.getProp<normal>(a)[0] = -1.0 * vd.getProp<normal>(a)[0];
                //vd.getProp<normal>(a)[1] = -1.0 * vd.getProp<normal>(a)[1];
            }
        }
        ++part;
    }
}


template <typename CellList> inline void density_summation(particles & vd, CellList & NN)
{
	vd.template ghost_get<rho>();
	auto part = vd.getDomainIterator();
	vd.updateCellList(NN);
	while (part.isNext())
	{
		auto a = part.get();
		Point<3, double> xa = vd.getPos(a);
		double rho_a = 0.0;
		double m = (return_sign(vd.getProp<sdf>(a)) < 0)?m_0:m_1;
		
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));
		
		while (Np.isNext() == true)
		{
			auto b = Np.get();
			Point<3, double> xb = vd.getPos(b);
			Point<3, double> dr = xa - xb;
			double r2 = norm2(dr);
			
			if (r2 <= 4.0*H*H)
			{
				double r = std::sqrt(r2);
				rho_a += Wab(r);
			}
			++Np;
		}
		vd.getProp<rho>(a) = m*rho_a;
		++part;
	}
}

inline void EqOfState(particles & vd)
{
	auto it = vd.getDomainIterator();
	while (it.isNext())
	{
		auto a = it.get();
		double ref_dens = 0.0;
		if (return_sign(vd.getProp<sdf>(a)) < 0) ref_dens = rho_0;
		else ref_dens = rho_1;
		double dens =  vd.getProp<rho>(a);
		vd.getProp<p>(a) = (c*c*ref_dens/gamma_eos)*(std::pow(dens/ref_dens, gamma_eos) - 1) + p_background;

		++it;
	}
}

template<typename CellList> inline void calc_forces(particles & vd, CellList & NN)
{
	vd.template ghost_get<rho, p, vel>();
	vd.updateCellList(NN);
	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	double max_p = 0.0;
	double max_eta = 0.0;
	double avg_alpha = 0.0;
	double maxvel = 0.0;
	int numparticles = 0;
	int numsurfparticles = 0;

	// For each particle ...
	while (part.isNext())
	{
		auto a = part.get();

		vd.getProp<force>(a)[0] = 0.0;
		vd.getProp<force>(a)[1] = 0.0;
		vd.getProp<force>(a)[2] = 0.0;

		Point<3, double> p_force;
		Point<3, double> eta_force;
		for(int k=0;k<3;k++) p_force[k] = 0.0;
		for(int k=0;k<3;k++) eta_force[k] = 0.0;

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

		// Get the Velocity of the particle
		Point<3,double> va = vd.getProp<vel>(a);

		// Get an iterator over the neighborhood particles of p
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

		// For each neighborhood particle
		while (Np.isNext() == true)
		{
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

				//double factor_p = - (Va*Va + Vb*Vb)*(rhoa*Pb + rhob*Pa)/(rhoa + rhob);
				double factor_p = -(Pa/(rhoa*rhoa) + Pb/(rhob*rhob))*massb;
				//double factor_p = - rhoa*massa*(Pa/(rhoa*rhoa) + Pb/(rhob*rhob))*massb;
				double factor_visc = (2*etaa*etab/(etaa + etab))*(Va*Va + Vb*Vb)*dwdrab/r;
				//if (t > 0.0002) factor_p = 0.0;

				p_force[0] += factor_p * DW.get(0);
				p_force[1] += factor_p * DW.get(1);
				p_force[2] += factor_p * DW.get(2);
				eta_force[0] += factor_visc*v_rel[0]/massa;
				eta_force[1] += factor_visc*v_rel[1]/massa;
				eta_force[2] += factor_visc*v_rel[2]/massa;

				//std::cout<<"W(0) = "<<Wab(0.0)<<" W(1.3/32) = "<<Wab(H)<<std::endl;
				//std::cout<<"Contribution by grad p: "<<factor_p*DW.get(2)<<", contribution by viscosity: "<<factor_visc*v_rel[2]<<", grad W magnitude: "<<DW.norm()<<", dwdr: "<<dwdrab<<std::endl;
			}

			++Np;
		}

		vd.getProp<force>(a)[0] += p_force[0] + eta_force[0];
		vd.getProp<force>(a)[1] += p_force[1] + eta_force[1];
		vd.getProp<force>(a)[2] += p_force[2] + eta_force[2];

		if (std::abs(vd.getProp<sdf>(a)) < (2.0*H))
		{
			double smoothing_factor = Wab(std::abs(vd.getProp<sdf>(a)))/(rhoa*Wab(0.0)); //latter is part of the equations, not the smoothing factor itself, but for convinience its in there
			//if (t < 0.1) smoothing_factor = std::sin(2*M_PI*t*2.5)*smoothing_factor;
			double sign_corr = 1.0;
			if (return_sign(vd.getProp<sdf>(a) > 0)) sign_corr = -1.0;
			const double stf0 = smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[0]*sign_corr;
			const double stf1 = smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[1]*sign_corr;
			const double stf2 = smoothing_factor*alpha*vd.getProp<curvature>(a)*vd.getProp<normal>(a)[2]*sign_corr;

			vd.getProp<force>(a)[0] += stf0;
			vd.getProp<force>(a)[1] += stf1;
			vd.getProp<force>(a)[2] += stf2;

			avg_alpha += vd.getProp<curvature>(a);
			numsurfparticles++;
			//if ((a.getKey() == 6342)||(a.getKey() == 5285)) std::cout<<"Particle "<<a.getKey()<<" Position: "<<vd.getPos(a)[0]<<", "<<vd.getPos(a)[1]<<", "<<vd.getPos(a)[2]<<", surface tension force: "<<stf0<<", "<<stf1<<", "<<stf2<<std::endl;
		}
		if (va.norm() > maxvel) maxvel = va.norm();
		if (p_force.norm() > max_p) max_p = p_force.norm();
		if (eta_force.norm() > max_eta) max_eta = eta_force.norm();
		++numparticles;
		++part;
	}
	avg_alpha = avg_alpha/numsurfparticles;
	if (corr) std::cout<<"Time step: "<<t<<", Max p force: "<<max_p<<", Max eta force: "<<max_eta<<", Max vel: "<<maxvel<<", Average curvature: "<<avg_alpha<<", number of particles: "<<numparticles<<std::endl;
	//if (corr) std::cout<<avg_alpha<<std::endl;
	if (numparticles == 0) throw std::invalid_argument("no particles");
}

void pred_corr_int(particles & vd, double dt, int & corr)
{
	// particle iterator
	auto part = vd.getDomainIterator();

	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		if (!corr)
		{
			// store the values at time n:
			vd.getProp<pos_prev>(a) = vd.getPos(a);
			vd.getProp<vel_prev>(a) = vd.getProp<vel>(a);
			// calculate intermediate values
			Point<3, double> dx = 0.5*dt*vd.getProp<vel>(a);
			vd.getPos(a)[0] += dx[0];
			vd.getPos(a)[1] += dx[1];
			vd.getPos(a)[2] += dx[2];
			vd.getProp<vel>(a) = vd.getProp<vel>(a) + 0.5*dt*vd.getProp<force>(a);

		}
		else
		{
		    //auto &v_cl=create_vcluster();
		    //std::cout<<v_cl.rank()<<": PosPrev:("<<vd.getProp<pos_prev>(a)[0]<<":"<<vd.getProp<pos_prev>(a)[1]<<":"<<vd.getProp<pos_prev>(a)[2]<<")"<<std::endl;
            //std::cout<<v_cl.rank()<<": VelPrev:("<<vd.getProp<vel_prev>(a)[0]<<":"<<vd.getProp<vel_prev>(a)[1]<<":"<<vd.getProp<vel_prev>(a)[2]<<")"<<std::endl;
            //std::cout<<v_cl.rank()<<": Force:("<<vd.getProp<force>(a)[0]<<":"<<vd.getProp<force>(a)[1]<<":"<<vd.getProp<force>(a)[2]<<")"<<std::endl;
			// correct the accelerations and velocities
			Point<3, double> x = vd.getProp<pos_prev>(a) + dt*(vd.getProp<vel_prev>(a) + 0.5*dt*vd.getProp<force>(a));
			vd.getPos(a)[0] = x[0];
			vd.getPos(a)[1] = x[1];
			vd.getPos(a)[2] = x[2];
			vd.getProp<vel>(a) = vd.getProp<vel_prev>(a) + dt*vd.getProp<force>(a);

		}
		++part;
	}
	corr = !corr;

}

int main(int argc, char* argv[])
{
	openfpm_init(&argc, &argv);

	Box<3, double> domain({-l_1/2.0, -l_1/2.0, -l_1/2.0}, {l_1/2.0, l_1/2.0, l_1/2.0});
	size_t sz[3] = {(size_t)(l_1/dp),(size_t)(l_1/dp),(size_t)(l_1/dp)};

	size_t bc[3] = {PERIODIC, PERIODIC, PERIODIC};
	Ghost<3, double> g(band_width);

	particles vd(0, domain, bc, g, DEC_GRAN(512));
	//particles_surface vd_s(vd.getDecomposition(), 0);

	openfpm::vector<std::string> names({"sdf","normal","curvature","surf_flag","rho","p","f","vel","vel_prev"});
	vd.setPropNames(names);

	auto particle_it = vd.getGridIterator(sz);

	while (particle_it.isNext())
	{
		double x = -l_1/2.0 + dp/2.0 + particle_it.get().get(0)*dp;
		double y = -l_1/2.0 + dp/2.0 + particle_it.get().get(1)*dp;
		double z = -l_1/2.0 + dp/2.0 + particle_it.get().get(2)*dp;

		vd.add();
		vd.getLastPos()[0] = x;
		vd.getLastPos()[1] = y;
		vd.getLastPos()[2] = z;
		vd.template getLastProp<surf_flag>() = 0;

		double dist = getDistancePointCube(abs(x), abs(y), abs(z));
		if ((abs(x) < l_0/2.0) && (abs(y) < l_0/2.0) && (abs(z) < l_0/2.0)) dist = -1.0*dist; // we're inside then

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
	vd.ghost_get<rho>();
	auto NN = vd.getCellList(2*H);
	density_summation(vd, NN);

	EqOfState(vd);

	vd.write("droplet_init_before_redistancing");

	Redist_options rdistoptions;
	rdistoptions.H = dp;
	rdistoptions.r_cutoff_factor = 2.4;
	rdistoptions.sampling_radius = 0.75*band_width;
	rdistoptions.tolerance = 1e-5;
	rdistoptions.compute_normals = 1;
	rdistoptions.compute_curvatures = 1;

	particle_cp_redistancing <particles, taylor4>pcprdist_init(vd, rdistoptions);

	//pcprdist_init.run_redistancing();

	vd.write("droplet_init");

	rdistoptions.tolerance = 1e-12;
	size_t write = 0;
	size_t rdist = 0;
	size_t rdist_small = 0;
	size_t it = 0;
	size_t it_reb = 0;
	t = 0.0;
	corr = 0;

	while (t <= t_end)
	{
		Vcluster<> & v_cl = create_vcluster();
		timer it_time;
		vd.map();
		vd.ghost_get<rho>();
		vd.updateCellList(NN);
		// Calculate pressure from the density
		density_summation(vd, NN);
		EqOfState(vd);

		if (rdist < t*1000)
		{
			if (t<0.05)	rdistoptions.write_sdf = 0;
			else rdistoptions.write_sdf = 1;
			//std::cout.setstate(std::ios_base::failbit);
			particle_cp_redistancing <particles, taylor4>pcprdist(vd, rdistoptions);
			pcprdist.run_redistancing();
			//std::cout.clear();
			rdist++;
		}

		// Calc forces
		calc_forces(vd, NN);

		// Predictor step
		it++;
		pred_corr_int(vd, dt, corr);
		vd.map();

		density_summation(vd, NN);
		EqOfState(vd);
		calc_forces(vd, NN);

		// corrector step
		pred_corr_int(vd, dt, corr);

		t += dt;


		if (write < t*100)
		{
			vd.deleteGhost();

			vd.write_frame("Geometry_droplet_v9/Geometry_droplet_rdist_v9",write);
			write++;

			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "TIME: " << t << std::endl;}

		}

	}
	openfpm_finalize();
}




