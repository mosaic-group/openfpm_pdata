/*!
 * \page Vector_7_sph_dlb_opt2 Vector 7 SPH Dam break simulation with Dynamic load balacing (Many core version)
 *
 *
 * [TOC]
 *
 *
 * # SPH with Dynamic load Balancing # {#SPH_dlb}
 *
 * This is just a rework of the SPH Dam break simulation optimized to get better performance. In this case we do not use the CRS scheme
 * But the nomeal symmetric scheme. This is better in case of many core. We do not provide further explanation because the optimization
 * operated here is the same as 
 *
 * \see \ref Vector_5_md_vl_sym
 *
 */

//#define SE_CLASS1
//#define STOP_ON_ERROR

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"
#include "util/stat/common_statistics.hpp"

// A constant to indicate boundary particles
#define BOUNDARY 1

// A constant to indicate fluid particles
#define FLUID 0

// initial spacing between particles dp in the formulas
const double dp = 0.00425;
// Maximum height of the fluid water
// is going to be calculated and filled later on
double h_swl = 0.0;

// c_s in the formulas (constant used to calculate the sound speed)
const double coeff_sound = 20.0;

// gamma in the formulas
const double gamma_ = 7.0;

// sqrt(3.0*dp*dp) support of the kernel
const double H = 0.00736121593217;

const double FourH2 = 4.0 * H*H;

// Eta in the formulas
const double Eta2 = 0.01 * H*H;

// alpha in the formula
const double visco = 0.1;

// cbar in the formula (calculated later)
double cbar = 0.0;

// Mass of the fluid particles
const double MassFluid = 0.0000767656;

// Mass of the boundary particles
const double MassBound = 0.0000767656;

// End simulation time
#ifndef TEST_RUN
const double t_end = 1.5;
#else
const double t_end = 0.004;
#endif

// Gravity acceleration
const double gravity = 9.81;

// Reference densitu 1000Kg/m^3
const double rho_zero = 1000.0;

// Filled later require h_swl, it is b in the formulas
double B = 0.0;

// Constant used to define time integration
const double CFLnumber = 0.2;

// Minimum T
const double DtMin = 0.00001;

// Minimum Rho allowed
const double RhoMin = 700.0;

// Maximum Rho allowed
const double RhoMax = 1300.0;

// Filled in initialization
double max_fluid_height = 0.0;

// Properties

// FLUID or BOUNDARY
const size_t type = 0;

// FLUID or BOUNDARY
const size_t nn_num = 1;

// Density
const int rho = 2;

// Density at step n-1
const int rho_prev = 3;

// Pressure
const int Pressure = 4;

// Delta rho calculated in the force calculation
const int drho = 5;

// calculated force
const int force = 6;

// velocity
const int velocity = 7;

// velocity at previous step
const int velocity_prev = 8;

/*! \cond [sim parameters] \endcond */

/*! \cond [vector_dist_def] \endcond */

// Type of the vector containing particles
typedef vector_dist<3,double,aggregate<int, int,double,  double,    double,     double,     double[3], double[3], double[3]> > particles;
//                                       |   |    |        |          |            |            |         |            |
//                                       |   |    |        |          |            |            |         |            |
//                                     type  |  density   density    Pressure    delta       force     velocity    velocity
//                                           |           at n-1                 density                           at n - 1
//                                           |
//									Number of neighborhood

//! Model for Dynamic load balancing
struct ModelCustom
{
	template<typename Decomposition, typename vector> inline void addComputation(Decomposition & dec, vector & vd, size_t v, size_t p)
	{
                if (vd.template getProp<type>(p) == FLUID )
                {dec.addComputationCost(v,4);}
                else
                {dec.addComputationCost(v,3);}
	}

	template<typename Decomposition> inline void applyModel(Decomposition & dec, size_t v)
	{
		dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
	}

    float distributionTol()
	{
		return 1.01;
	}
};

//! Second model for dynamic load balancing
struct ModelCustom2
{
	template<typename Decomposition, typename vector> inline void addComputation(Decomposition & dec, vector & vd, size_t v, size_t p)
	{
		dec.addComputationCost(v,vd.template getProp<nn_num>(p) + 4);
	}

	template<typename Decomposition> inline void applyModel(Decomposition & dec, size_t v)
	{
	}

	float distributionTol()
	{
			return 1.01;
	}
};

inline void EqState(particles & vd)
{
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto a = it.get();

		double rho_a = vd.template getProp<rho>(a);
		double rho_frac = rho_a / rho_zero;

		vd.template getProp<Pressure>(a) = B*( rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac - 1.0);

		++it;
	}
}

const double a2 = 1.0/M_PI/H/H/H;

inline double Wab(double r)
{
	r /= H;

	if (r < 1.0)
	{return (1.0 - 3.0/2.0*r*r + 3.0/4.0*r*r*r)*a2;}
	else if (r < 2.0)
	{return (1.0/4.0*(2.0 - r*r)*(2.0 - r*r)*(2.0 - r*r))*a2;}
	else
	{return 0.0;}
}

const double c1 = -3.0/M_PI/H/H/H/H;
const double d1 = 9.0/4.0/M_PI/H/H/H/H;
const double c2 = -3.0/4.0/M_PI/H/H/H/H;
const double a2_4 = 0.25*a2;
// Filled later
double W_dap = 0.0;

inline void DWab(Point<3,double> & dx, Point<3,double> & DW, double r)
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
}

inline double Tensile(double r, double rhoa, double rhob, double prs1, double prs2)
{
	const double qq=r/H;
	//-Cubic Spline kernel
	double wab;
	if(r>H)
	{
		double wqq1=2.0f-qq;
		double wqq2=wqq1*wqq1;

		wab=a2_4*(wqq2*wqq1);
	}
	else
	{
	    double wqq2=qq*qq;
	    double wqq3=wqq2*qq;

	    wab=a2*(1.0f-1.5f*wqq2+0.75f*wqq3);
	}

	//-Tensile correction.
	double fab=wab*W_dap;
	fab*=fab; fab*=fab; //fab=fab^4
	const double tensilp1=(prs1/(rhoa*rhoa))*(prs1>0? 0.01: -0.2);
	const double tensilp2=(prs2/(rhob*rhob))*(prs2>0? 0.01: -0.2);

	return (fab*(tensilp1+tensilp2));
}


inline double Pi(const Point<3,double> & dr, double rr2, Point<3,double> & dv, double rhoa, double rhob, double massb, double & visc)
{
	const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
	const double dot_rr2 = dot/(rr2+Eta2);
	visc=std::max(dot_rr2,visc);

	if(dot < 0)
	{
		const float amubar=H*dot_rr2;
		const float robar=(rhoa+rhob)*0.5f;
		const float pi_visc=(-visco*cbar*amubar/robar);

		return pi_visc;
    }
	else
		return 0.0;
}


template<typename VerletList> inline void calc_forces(particles & vd, VerletList & NN, double & max_visc, timer & tf)
{
	/*! \cond [reset_particles] \endcond */

	// Reset the ghost
    auto itg = vd.getDomainIterator();
    while (itg.isNext())
    {
        auto p = itg.get();
        // Reset force

		// Reset the force counter (- gravity on zeta direction)
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = -gravity;
		vd.template getProp<drho>(p) = 0.0;


        ++itg;
    }

    /*! \cond [reset_particles] \endcond */

    /*! \cond [reset_particles2] \endcond */

    auto itg2 = vd.getGhostIterator();
    while (itg2.isNext())
    {
        auto p = itg2.get();
        // Reset force

		// Reset the force counter (- gravity on zeta direction)
		vd.template getProp<force>(p)[0] = 0.0;
		vd.template getProp<force>(p)[1] = 0.0;
		vd.template getProp<force>(p)[2] = 0.0;
		vd.template getProp<drho>(p) = 0.0;

        ++itg2;
    }

    /*! \cond [reset_particles2] \endcond */

	tf.start();

    // Get an iterator over particles
    auto part = vd.getDomainIterator();

	double visc = 0;

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// Get the position xp of the particle
		Point<3,double> xa = vd.getPos(a);

		// Type of the particle
		size_t typea = vd.getProp<type>(a);

		// Take the mass of the particle dependently if it is FLUID or BOUNDARY
		double massa = (vd.getProp<type>(a) == FLUID)?MassFluid:MassBound;

		// Get the density of the of the particle a
		double rhoa = vd.getProp<rho>(a);

		// Get the pressure of the particle a
		double Pa = vd.getProp<Pressure>(a);

		// Get the Velocity of the particle a
		Point<3,double> va = vd.getProp<velocity>(a);

        // Get an iterator over the neighborhood particles of p
        auto Np = NN.template getNNIterator<NO_CHECK>(a.getKey());

        size_t nn = 0;

        // For each neighborhood particle
        while (Np.isNext() == true)
        {
                // ... q
                auto b = Np.get();

                // Get the position xp of the particle
                Point<3,double> xb = vd.getPos(b);

                // Get the distance between p and q
                Point<3,double> dr = xa - xb;
                // take the norm of this vector
                float r2 = norm2(dr);

                // if they interact
                if (r2 < FourH2 && r2 > 1e-18)
                {
                        double r = sqrt(r2);

                        unsigned int typeb = vd.getProp<type>(b);

                        double massb = (typeb == FLUID)?MassFluid:MassBound;
                        Point<3,double> vb = vd.getProp<velocity>(b);
                        double Pb = vd.getProp<Pressure>(b);
                        double rhob = vd.getProp<rho>(b);

                        Point<3,double> v_rel = va - vb;

                        Point<3,double> DW;
                        DWab(dr,DW,r);

                        //! \cond [symmetry2] \endcond

                        double factor = - ((Pa + Pb) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb) + Pi(dr,r2,v_rel,rhoa,rhob,massb,visc));

                        // Bound - Bound does not produce any change
                        factor = (typea == BOUNDARY && typeb == BOUNDARY)?0.0f:factor;

                        vd.getProp<force>(a)[0] += massb * factor * DW.get(0);
                        vd.getProp<force>(a)[1] += massb * factor * DW.get(1);
                        vd.getProp<force>(a)[2] += massb * factor * DW.get(2);

                        vd.getProp<force>(b)[0] -= massa * factor * DW.get(0);
                        vd.getProp<force>(b)[1] -= massa * factor * DW.get(1);
                        vd.getProp<force>(b)[2] -= massa * factor * DW.get(2);

                        double scal = (v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2));

                        // Bound - Bound does not produce any change
                        scal = (typea == BOUNDARY && typeb == BOUNDARY)?0.0f:scal;

                        vd.getProp<drho>(a) += massb*scal;
                        vd.getProp<drho>(b) += massa*scal;

                        //! \cond [symmetry2] \endcond
                }

                nn++;
                ++Np;
        }

        // Number of particles here
        vd.getProp<nn_num>(a) = nn;

		++part;
	}

	tf.stop();

	//! \cond [ghost_put] \endcond

	vd.template ghost_put<add_,drho,force>();

	//! \cond [ghost_put] \endcond
}

void max_acceleration_and_velocity(particles & vd, double & max_acc, double & max_vel)
{
	// Calculate the maximum acceleration
	auto part = vd.getDomainIterator();

	while (part.isNext())
	{
		auto a = part.get();

		Point<3,double> acc(vd.getProp<force>(a));
		double acc2 = norm2(acc);

		Point<3,double> vel(vd.getProp<velocity>(a));
		double vel2 = norm2(vel);

		if (vel2 >= max_vel)
			max_vel = vel2;

		if (acc2 >= max_acc)
			max_acc = acc2;

		++part;
	}
	max_acc = sqrt(max_acc);
	max_vel = sqrt(max_vel);

	Vcluster<> & v_cl = create_vcluster();
	v_cl.max(max_acc);
	v_cl.max(max_vel);
	v_cl.execute();
}

double calc_deltaT(particles & vd, double ViscDtMax)
{
	double Maxacc = 0.0;
	double Maxvel = 0.0;
	max_acceleration_and_velocity(vd,Maxacc,Maxvel);

	//-dt1 depends on force per unit mass.
	const double dt_f = (Maxacc)?sqrt(H/Maxacc):std::numeric_limits<int>::max();

	//-dt2 combines the Courant and the viscous time-step controls.
	const double dt_cv = H/(std::max(cbar,Maxvel*10.) + H*ViscDtMax);

	//-dt new value of time step.
	double dt=double(CFLnumber)*std::min(dt_f,dt_cv);
	if(dt<double(DtMin))
	{dt=double(DtMin);}

	return dt;
}


openfpm::vector<size_t> to_remove;

size_t cnt = 0;

/*! \cond [verlet_new_arg] \endcond */
void verlet_int(particles & vd, double dt, double & max_disp)
/*! \cond [verlet_new_arg] \endcond */
{
	// list of the particle to remove
	to_remove.clear();

	// particle iterator
	auto part = vd.getDomainIterator();

	double dt205 = dt*dt*0.5;
	double dt2 = dt*2.0;

	/*! \cond [reset_max_disp] \endcond */
	max_disp = 0;
	/*! \cond [reset_max_disp] \endcond */

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// if the particle is boundary
		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			// Update rho
			double rhop = vd.template getProp<rho>(a);

			// Update only the density
	    	vd.template getProp<velocity>(a)[0] = 0.0;
	    	vd.template getProp<velocity>(a)[1] = 0.0;
	    	vd.template getProp<velocity>(a)[2] = 0.0;
            double rhonew = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);
            vd.template getProp<rho>(a) = (rhonew < rho_zero)?rho_zero:rhonew;

		    vd.template getProp<rho_prev>(a) = rhop;

			++part;
			continue;
		}

		//-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
		double dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
	    double dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
	    double dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;

	    vd.getPos(a)[0] += dx;
	    vd.getPos(a)[1] += dy;
	    vd.getPos(a)[2] += dz;

	    /*! \cond [calc_max_disp] \endcond */
	    double d2 = dx*dx + dy*dy + dz*dz;

	    max_disp = (max_disp > d2)?max_disp:d2;
	    /*! \cond [calc_max_disp] \endcond */

	    double velX = vd.template getProp<velocity>(a)[0];
	    double velY = vd.template getProp<velocity>(a)[1];
	    double velZ = vd.template getProp<velocity>(a)[2];
	    double rhop = vd.template getProp<rho>(a);

    	vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<force>(a)[0]*dt2;
    	vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<force>(a)[1]*dt2;
    	vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<force>(a)[2]*dt2;
    	vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);

	    // Check if the particle go out of range in space and in density
	    if (vd.getPos(a)[0] <  0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
	        vd.getPos(a)[0] >  0.000263878+1.59947 || vd.getPos(a)[1] > 0.000263878+0.672972 || vd.getPos(a)[2] > 0.000263878+0.903944 ||
			vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
	    {
	                   to_remove.add(a.getKey());
	    }

	    vd.template getProp<velocity_prev>(a)[0] = velX;
	    vd.template getProp<velocity_prev>(a)[1] = velY;
	    vd.template getProp<velocity_prev>(a)[2] = velZ;
	    vd.template getProp<rho_prev>(a) = rhop;

		++part;
	}

	/*! \cond [max_across_proc] \endcond */

	Vcluster<> & v_cl = create_vcluster();
	v_cl.max(max_disp);
	v_cl.execute();

	max_disp = sqrt(max_disp);

	/*! \cond [max_across_proc] \endcond */

	// increment the iteration counter
	cnt++;
}

/*! \cond [euler_new_arg] \endcond */
void euler_int(particles & vd, double dt, double & max_disp)
/*! \cond [euler_new_arg] \endcond */
{
	// list of the particle to remove
	to_remove.clear();

	// particle iterator
	auto part = vd.getDomainIterator();

	double dt205 = dt*dt*0.5;
	double dt2 = dt*2.0;

	max_disp = 0;

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// if the particle is boundary
		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			// Update rho
			double rhop = vd.template getProp<rho>(a);

			// Update only the density
	    	vd.template getProp<velocity>(a)[0] = 0.0;
	    	vd.template getProp<velocity>(a)[1] = 0.0;
	    	vd.template getProp<velocity>(a)[2] = 0.0;
            double rhonew = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);
            vd.template getProp<rho>(a) = (rhonew < rho_zero)?rho_zero:rhonew;

		    vd.template getProp<rho_prev>(a) = rhop;

			++part;
			continue;
		}

		//-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
		double dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
	    double dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
	    double dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;

	    vd.getPos(a)[0] += dx;
	    vd.getPos(a)[1] += dy;
	    vd.getPos(a)[2] += dz;

	    double d2 = dx*dx + dy*dy + dz*dz;

	    max_disp = (max_disp > d2)?max_disp:d2;

	    double velX = vd.template getProp<velocity>(a)[0];
	    double velY = vd.template getProp<velocity>(a)[1];
	    double velZ = vd.template getProp<velocity>(a)[2];
	    double rhop = vd.template getProp<rho>(a);

    	vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + vd.template getProp<force>(a)[0]*dt;
    	vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + vd.template getProp<force>(a)[1]*dt;
	   	vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + vd.template getProp<force>(a)[2]*dt;
	   	vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);

	    // Check if the particle go out of range in space and in density
	    if (vd.getPos(a)[0] <  0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
	        vd.getPos(a)[0] >  0.000263878+1.59947 || vd.getPos(a)[1] > 0.000263878+0.672972 || vd.getPos(a)[2] > 0.000263878+0.903944 ||
			vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
	    {
	                   to_remove.add(a.getKey());
	    }

	    vd.template getProp<velocity_prev>(a)[0] = velX;
	    vd.template getProp<velocity_prev>(a)[1] = velY;
	    vd.template getProp<velocity_prev>(a)[2] = velZ;
	    vd.template getProp<rho_prev>(a) = rhop;

		++part;
	}

	Vcluster<> & v_cl = create_vcluster();
	v_cl.max(max_disp);
	v_cl.execute();

	max_disp = sqrt(max_disp);

	// increment the iteration counter
	cnt++;
}

int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<3,double> domain({-0.05,-0.05,-0.05},{1.7010,0.7065,0.5025});
	size_t sz[3] = {413,179,133};

	// Fill W_dap
	W_dap = 1.0/Wab(H/1.5);

	// Here we define the boundary conditions of our problem
    size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

    /*! \cond [skin_calc] \endcond */

    double skin = 0.25 * 2*H;
    double r_gskin = 2*H + skin;

	// extended boundary around the domain, and the processor domain
    // by the support of the cubic kernel
	Ghost<3,double> g(r_gskin);

	/*! \cond [important_option] \endcond */

	particles vd(0,domain,bc,g,DEC_GRAN(512));

	/*! \cond [important_option] \endcond */

	// You can ignore all these dp/2.0 is a trick to reach the same initialization
	// of Dual-SPH that use a different criteria to draw particles
	Box<3,double> fluid_box({dp/2.0,dp/2.0,dp/2.0},{0.4+dp/2.0,0.67-dp/2.0,0.3+dp/2.0});

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd,sz,domain,fluid_box);

	// here we fill some of the constants needed by the simulation
	max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
	h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
	B = (coeff_sound)*(coeff_sound)*gravity*h_swl*rho_zero / gamma_;
	cbar = coeff_sound * sqrt(gravity * h_swl);

	// for each particle inside the fluid box ...
        while (fluid_it.isNext())
	{
		// ... add a particle ...
		vd.add();

		// ... and set it position ...
		vd.getLastPos()[0] = fluid_it.get().get(0);
		vd.getLastPos()[1] = fluid_it.get().get(1);
		vd.getLastPos()[2] = fluid_it.get().get(2);

		// and its type.
		vd.template getLastProp<type>() = FLUID;

		// We also initialize the density of the particle and the hydro-static pressure given by
		//
		// rho_zero*g*h = P
		//
		// rho_p = (P/B + 1)^(1/Gamma) * rho_zero
		//

		vd.template getLastProp<Pressure>() = rho_zero * gravity *  (max_fluid_height - fluid_it.get().get(2));

		vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0/gamma_) * rho_zero;
		vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		// next fluid particle
		++fluid_it;
        }

	// Recipient
	Box<3,double> recipient1({0.0,0.0,0.0},{1.6+dp/2.0,0.67+dp/2.0,0.4+dp/2.0});
	Box<3,double> recipient2({dp,dp,dp},{1.6-dp/2.0,0.67-dp/2.0,0.4+dp/2.0});

	Box<3,double> obstacle1({0.9,0.24-dp/2.0,0.0},{1.02+dp/2.0,0.36,0.45+dp/2.0});
	Box<3,double> obstacle2({0.9+dp,0.24+dp/2.0,0.0},{1.02-dp/2.0,0.36-dp,0.45-dp/2.0});
	Box<3,double> obstacle3({0.9+dp,0.24,0.0},{1.02,0.36,0.45});

	openfpm::vector<Box<3,double>> holes;
	holes.add(recipient2);
	holes.add(obstacle1);
	auto bound_box = DrawParticles::DrawSkin(vd,sz,domain,holes,recipient1);

	while (bound_box.isNext())
	{
		vd.add();

		vd.getLastPos()[0] = bound_box.get().get(0);
		vd.getLastPos()[1] = bound_box.get().get(1);
		vd.getLastPos()[2] = bound_box.get().get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = rho_zero;
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		++bound_box;
	}

	auto obstacle_box = DrawParticles::DrawSkin(vd,sz,domain,obstacle2,obstacle1);

	while (obstacle_box.isNext())
	{
		vd.add();

		vd.getLastPos()[0] = obstacle_box.get().get(0);
		vd.getLastPos()[1] = obstacle_box.get().get(1);
		vd.getLastPos()[2] = obstacle_box.get().get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = rho_zero;
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		++obstacle_box;
	}

	vd.map();

	//////// Number of particles
	//
	//
	//

	size_t tot_part = vd.size_local();

	auto & v_cl = create_vcluster();
	v_cl.sum(tot_part);
	v_cl.execute();
	std::cout << "SUM: " << tot_part << std::endl;

	// Now that we fill the vector with particles
	ModelCustom md;

	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map();

	vd.ghost_get<type,rho,Pressure,velocity>();

	/*! \cond [get_verlet_sym] \endcond */
	auto NN = vd.getVerletSym(r_gskin);
	/*! \cond [get_verlet_sym] \endcond */

	openfpm::vector<double> time_forces;
	openfpm::vector<double> time_comm;
	openfpm::vector<double> time_steps;

	size_t write = 0;
	size_t it = 0;
	size_t it_reb = 0;
	double t = 0.0;
	double tot_disp = 0.0;
	double max_disp;
	while (t <= t_end)
	{
		Vcluster<> & v_cl = create_vcluster();
		timer it_time;

		it_time.start();

		/*! \cond [update_verlet] \endcond */

		it_reb++;
		if (2*tot_disp >= skin)
		{
			vd.remove(to_remove);

			vd.map();

			if (it_reb > 200)
			{
				ModelCustom2 md;
				vd.addComputationCosts(md);
				vd.getDecomposition().redecompose(200);

				vd.map();

				it_reb = 0;

				if (v_cl.getProcessUnitID() == 0)
					std::cout << "REBALANCED " << std::endl;

			}

			// Calculate pressure from the density
			EqState(vd);

			vd.ghost_get<type,rho,Pressure,velocity>();

			/*! \cond [update_verlet_sym] \endcond */
			vd.updateVerlet(NN,r_gskin,VL_SYMMETRIC);
			/*! \cond [update_verlet_sym] \endcond */

			tot_disp = 0.0;

			if (v_cl.getProcessUnitID() == 0)
				std::cout << "RECONSTRUCT Verlet " << std::endl;
		}
		else
		{
			// Calculate pressure from the density
			EqState(vd);

			timer tc;
			tc.start();

			vd.ghost_get<type,rho,Pressure,velocity>(SKIP_LABELLING);

			tc.stop();
			time_comm.add(tc.getwct());
		}

		/*! \cond [update_verlet] \endcond */

		double max_visc = 0.0;

		timer tf;

		// Calc forces
		calc_forces(vd,NN,max_visc,tf);

		// Get the maximum viscosity term across processors
		v_cl.max(max_visc);
		v_cl.execute();

		// Calculate delta t integration
		double dt = calc_deltaT(vd,max_visc);

		/*! \cond [pass_ver_eu] \endcond */

		// VerletStep or euler step
		it++;
		if (it < 40)
			verlet_int(vd,dt,max_disp);
		else
		{
			euler_int(vd,dt,max_disp);
			it = 0;
		}

		tot_disp += max_disp;

		/*! \cond [pass_ver_eu] \endcond */

		t += dt;

		it_time.stop();

		if (write < t*100)
		{
			if (v_cl.getProcessUnitID() == 0)
				std::cout << "TIME: " << t << " " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "  TOT disp: " << tot_disp << "    " << cnt << std::endl;

		}
		else
		{
			if (v_cl.getProcessUnitID() == 0)
				std::cout << "TIME: " << t << "  " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "  TOT disp: " << tot_disp << "    " << cnt << std::endl;
		}
		time_steps.add(it_time.getwct());
		time_forces.add(tf.getwct());
	}

	double mean_ts, dev_ts;
	double mean_tf, dev_tf;
	double mean_comm, dev_comm;

	standard_deviation(time_steps,mean_ts,dev_ts);
	standard_deviation(time_forces,mean_tf,dev_tf);
	standard_deviation(time_comm,mean_comm,dev_comm);

	v_cl.sum(mean_ts);
	v_cl.sum(mean_tf);
	v_cl.sum(mean_comm);

	v_cl.sum(dev_ts);
	v_cl.sum(dev_tf);
	v_cl.sum(dev_comm);

	mean_ts /= v_cl.size();
	mean_tf /= v_cl.size();
	mean_comm /= v_cl.size();

	dev_ts /= v_cl.size();
	dev_tf /= v_cl.size();
	dev_comm /= v_cl.size();

	std::cout << mean_ts << " " << dev_ts << std::endl;
	std::cout << mean_tf << " " << dev_tf << std::endl;
	std::cout << mean_comm << " " << dev_comm << std::endl;

	openfpm_finalize();
}
 
