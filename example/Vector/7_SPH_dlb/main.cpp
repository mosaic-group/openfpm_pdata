/*!
 * \page Vector_7_sph_dlb Vector 7 SPH Dam break  simulation with Dynamic load balacing
 *
 *
 * [TOC]
 *
 *
 * # SPH with Dynamic load Balancing # {#SPH_dlb}
 *
 *
 * This example show the classical SPH Dam break simulation with Dynamic load balancing
 *
 * ## inclusion ## {#e0_v_inclusion}
 *
 * In order to use distributed vectors in our code we have to include the file Vector/vector_dist.hpp
 *
 * \snippet Vector/7_SPH_dlb/main.cpp inclusion
 *
 */

//! \cond [inclusion] \endcond
#include "Vector/vector_dist.hpp"
#include <math.h>
//! \cond [inclusion] \endcond

#include "Draw/DrawParticles.hpp"

#define BOUNDARY 0
#define FLUID 1

double lower_z = 2.0;

const double dp = 0.0085;
double h_swl = 0.0;
const double coeff_sound = 20.0;
const double gamma_ = 7.0;
// sqrt(3.0*dp*dp)
const double H = 0.0147224318643;
const double Eta2 = 0.01 * H*H;
const double visco = 0.1;
double cbar = 0.0;
const double MassFluid = 0.000614125;
const double MassBound = 0.000614125;
const double t_end = 1.0;
const double gravity = 9.81;
const double rho_zero = 1000.0;
double B = 0.0;
const double CFLnumber = 0.2;
const double DtMin = 0.00001;

// Filled in initialization
double max_fluid_height = 0.0;

const size_t type = 0;
const int rho = 1;
const int rho_prev = 2;
const int Pressure = 3;
const int drho = 4;
const int force = 5;
const int velocity = 6;
const int velocity_prev = 7;

typedef vector_dist<3,double,aggregate<size_t,double,double,double,double,double[3],double[3],double[3]>> particles;

inline void EqState(particles & vd)
{
//	double min = 100000000.0;
//	double max = 0.0;
//	double accum = 0.0;
//	size_t n = 0;

//	Vcluster & v_cl = create_vcluster();

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto a = it.get();

		double rho_a = vd.template getProp<rho>(a);
		double rho_frac = rho_a / rho_zero;

		vd.template getProp<Pressure>(a) = B*( rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac - 1.0);

		/// DEBUG

/*		if (vd.template getProp<Pressure>(a) < min)
			min = vd.template getProp<Pressure>(a);

		if (vd.template getProp<Pressure>(a) > max)
			max = vd.template getProp<Pressure>(a);

		if (vd.template getProp<Pressure>(a) > 2849.0)
			std::cout << "Particle: " << Point<3,double>(vd.getPos(a)).toString() << std::endl;

		accum += vd.template getProp<Pressure>(a);
		n++;*/

		++it;
	}

/*	v_cl.max(max);
	v_cl.min(min);
	v_cl.sum(accum);
	v_cl.sum(n);

	v_cl.execute();

	std::cout << "Max: " << max << " min: " << min << "  accum: " << accum/n << "     " << B << " n: " << n << std::endl;*/
}

const double c1 = -3.0/M_PI/H/H/H/H;
const double d1 = 9.0/4.0/M_PI/H/H/H/H;
const double c2 = -3.0/4.0/M_PI/H/H/H/H;
const double a2 = 1.0/M_PI/H/H/H;
const double a2_4 = 0.25*a2;
// Filled later
double W_dap = 0.0;

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


inline void DWab(Point<3,double> & dx, Point<3,double> & DW, double r, bool print)
{
	const double qq=r/H;

	if (qq < 1.0)
	{
		double qq2 = qq * qq;
		double fac = (c1*qq + d1*qq2)/r;

		DW.get(0) = fac*dx.get(0);
		DW.get(1) = fac*dx.get(1);
		DW.get(2) = fac*dx.get(2);
	}
	else if (qq < 2.0)
	{
		double wqq = (2.0 - qq);
		double fac = c2 * wqq * wqq / r;

		DW.get(0) = fac * dx.get(0);
		DW.get(1) = fac * dx.get(1);
		DW.get(2) = fac * dx.get(2);
	}
	else
	{
		DW.get(0) = 0.0;
		DW.get(1) = 0.0;
		DW.get(2) = 0.0;
	}
}

// Tensile correction
inline double Tensile(double r, double rhoa, double rhob, double prs1, double prs2, bool print)
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

//	if (print == true)
//		std::cout << "fab " << fab << " tensilp1: " << tensilp1 << " tensilp2: " << tensilp2 << std::endl;

	return (fab*(tensilp1+tensilp2));
}

inline double Pi(const Point<3,double> & dr, double rr2, Point<3,double> & dv, double rhoa, double rhob, double massb, double & visc, bool print)
{
	const double dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
	const double dot_rr2 = dot/(rr2+Eta2);
	visc=std::max(dot_rr2,visc);

	if(dot < 0)
	{
		const float amubar=H*dot_rr2;
		const float robar=(rhoa+rhob)*0.5f;
		const float pi_visc=(-visco*cbar*amubar/robar);

		if (print == true)
			std::cout << "   visco: " << visco << "  " << cbar << "  " << amubar << "   " << robar << "   ";

		return pi_visc;
    }
	else
		return 0.0;
}

template<typename CellList> inline double calc_forces(particles & vd, CellList & NN, double & max_visc)
{
	auto part = vd.getDomainIterator();

	Point<3,double> ForceMax({       0.0,      0.0,      0.0});
	Point<3,double> ForceMin({10000000.0,1000000.0,1000000,0});

	double visc = 0;
//	std::cout << "c1: " << c1 << "    c2: " << c2 << "    d1: " << d1 << std::endl;

	vd.updateCellList(NN);

	while (part.isNext())
	{
		auto a = part.get();

		// Get the position xp of the particle
		Point<3,double> xa = vd.getPos(a);

		if (vd.getProp<type>(a) != FLUID)
		{
			++part;
			continue;
		}

		double massa = (vd.getProp<type>(a) == FLUID)?MassFluid:MassBound;
		double rhoa = vd.getProp<rho>(a);
		double Pa = vd.getProp<Pressure>(a);
		Point<3,double> va = vd.getProp<velocity>(a);

		// Reset the force counter
		vd.template getProp<force>(a)[0] = 0.0;
		vd.template getProp<force>(a)[1] = 0.0;
		vd.template getProp<force>(a)[2] = -gravity;
		vd.template getProp<drho>(a) = 0.0;

		size_t cnt = 0;

//		std::cout << "---------------------" << std::endl;
//		std::cout << vd.getPos(a)[0] << "   " << vd.getPos(a)[1] << "   " << vd.getPos(a)[2] << std::endl;

		// Get an iterator over the neighborhood particles of p
		auto Np = NN.template getNNIterator<NO_CHECK>(NN.getCell(vd.getPos(a)));

		while (Np.isNext() == true)
		{
			// ... q
			auto b = Np.get();

			// Get the position xp of the particle
			Point<3,double> xb = vd.getPos(b);

			// if (p == q) skip this particle
			if (a.getKey() == b)	{++Np; continue;};

			double massb = (vd.getProp<type>(b) == FLUID)?MassFluid:MassBound;
			Point<3,double> vb = vd.getProp<velocity>(b);
			double Pb = vd.getProp<Pressure>(b);
			double rhob = vd.getProp<rho>(b);

			// Get the distance between p and q
			Point<3,double> dr = xa - xb;
			// take the norm of this vector
			double r2 = norm2(dr);

			if (r2 < 4.0*H*H)
			{
				double r = sqrt(r2);

				Point<3,double> v_rel = va - vb;

				Point<3,double> DW;
				DWab(dr,DW,r,false);

				double factor = - massb*((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb,false) + Pi(dr,r2,v_rel,rhoa,rhob,massb,visc,false));
/*
				if (create_vcluster().getProcessUnitID() == 0)
				{
					std::cout << "PARTICLE: " << dr.toString() << std::endl;
					std::cout << "Pressure: " << Pa << "  " << Pb << std::endl;
					std::cout << "Density: " << rhoa << "  " << rhob << std::endl;
					std::cout << "FACTOR: " << factor << std::endl;
					std::cout << "Tensile: " << Tensile(r,rhoa,rhob,Pa,Pb,false) << std::endl;
					std::cout << "DW: " << DW.get(0) << "     " << DW.get(1) << "     " << DW.get(2) << std::endl;
				}*/

				vd.getProp<force>(a)[0] += factor * DW.get(0);
				vd.getProp<force>(a)[1] += factor * DW.get(1);
				vd.getProp<force>(a)[2] += factor * DW.get(2);

				if (xa.get(0) > 0.0085 && xa.get(0) < 0.0105 &&
					xa.get(1) > 0.0085 && xa.get(1) < 0.0105 &&
					xa.get(2) > 0.0085 && xa.get(2) < 0.0105)
				{
					std::cout << "POSITION: " << xb.toString() << std::endl;
					std::cout << "FORCE: " << factor*DW.get(0) << "    " << factor*DW.get(1) << "      " << factor*DW.get(2) << std::endl;
					std::cout << "FORCE TERM: " << Tensile(r,rhoa,rhob,Pa,Pb,true) << "      " << Pi(dr,r2,v_rel,rhoa,rhob,massb,visc,true) << std::endl;
					cnt++;
				}

	            vd.getProp<drho>(a) += massb*(v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2));

/*	            if (vd.getProp<type>(a) == FLUID)
	            {
	            	std::cout << "DELTA RHO: " << massb*(v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2)) << "    " << v_rel.get(0) <<  "     " << v_rel.get(1) << "     "  << v_rel.get(2) <<  "  VISC: " << Pi(dr,r2,v_rel,rhoa,rhob,massb,visc,true) << std::endl;
	            }*/
			}

			++Np;
		}

		if (xa.get(0) > 0.0085 && xa.get(0) < 0.0105 &&
			xa.get(1) > 0.0085 && xa.get(1) < 0.0105 &&
			xa.get(2) > 0.0085 && xa.get(2) < 0.0105)
		{
			std::cout << "FORCE FINAL: " << vd.getProp<force>(a)[0] << "   " << vd.getProp<force>(a)[1] << "  " << vd.getProp<force>(a)[2] << "    " << cnt << std::endl;
		}

/*        if (vd.getProp<type>(a) == FLUID)
        {
        	std::cout << "DELTA DENSITY: " << vd.getProp<drho>(a) << std::endl;
        	std::cout << "+++++++++++++++++++++++++++++++++++" << std::endl;
        }*/

		++part;

		if (Point<3,double>(vd.getProp<force>(a)).norm() > ForceMax.norm()  )
		{
			ForceMax = Point<3,double>(vd.getProp<force>(a));
//			std::cout << "ForceMax: " << ForceMax.toString() << "   " << a.getKey() << std::endl;

			Point<3,double> p({0.01,0.0,0.0});
			Point<3,double> DW;

			DWab(p,DW,0.01,false);

//			std::cout << DW.get(0) << "   " << DW.get(1) << "      " << DW.get(2)  << std::endl;
		}

		if (Point<3,double>(vd.getProp<force>(a)).norm() < ForceMin.norm()  )
		{
			ForceMin = Point<3,double>(vd.getProp<force>(a));
		}
	}

	// Get the maximum viscosity term across processors
	Vcluster & v_cl = create_vcluster();
	v_cl.max(visc);
	v_cl.execute();
	max_visc = visc;

//	std::cout << "---------------------------------" << std::endl;
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
}

double dt_old = 0.0;

double calc_deltaT(particles & vd, double ViscDtMax)
{
	double Maxacc = 0.0;
	double Maxvel = 0.0;
	max_acceleration_and_velocity(vd,Maxacc,Maxvel);

	//-dt1 depends on force per unit mass.
	const double dt_f = (Maxacc)?sqrt(H/Maxacc):std::numeric_limits<int>::max();

	//-dt2 combines the Courant and the viscous time-step controls.
	const double dt_cv = H/(std::max(cbar,Maxvel*10.) + H*ViscDtMax);

//	std::cout << "dt_f   " << dt_f << "     dt_cv: " << dt_cv << "    Maxvel: " << Maxvel << std::endl;

	//-dt new value of time step.
	double dt=double(CFLnumber)*std::min(dt_f,dt_cv);
	if(dt<double(DtMin))
		dt=double(DtMin);

	if (dt_old != dt)
	{
		std::cout << "Dt changed to " << dt << "    MaxVel: " << Maxvel << "     Maxacc: " << Maxacc << "     lower_z: " << lower_z << std::endl;

		dt_old = dt;
	}

//	std::cout << "Returned dt: " << dt << std::endl;
	return dt;
}

openfpm::vector<size_t> to_remove;

bool particle_out = true;

void verlet_int(particles & vd, double dt)
{
	to_remove.clear();

	// Calculate the maximum acceleration
	auto part = vd.getDomainIterator();

	double dt205 = dt*dt*0.5;
	double dt2 = dt*2.0;

	while (part.isNext())
	{
		auto a = part.get();

		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			// Update only the density
			++part;
			continue;
		}

		//-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
		double dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
	    double dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
	    double dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;

//	    bool outrhop=(rhopnew<RhopOutMin||rhopnew>RhopOutMax);

//	    if (vd.getPos(a)[0] > 0.0085 && vd.getPos(a)[0] < 0.0105 &&
//	    	vd.getPos(a)[1] > 0.0085 && vd.getPos(a)[1] < 0.0105 &&
//			vd.getPos(a)[2] > 0.0085 && vd.getPos(a)[2] < 0.0105)
//	    {
//	    	std::cout << "DeltaX: " << dx << "    " << dy << "      " << dz << std::endl;
/*	    	std::cout << "FORCE: " << vd.template getProp<force>(a)[0] << "       "
	    			  << vd.template getProp<force>(a)[1] << "      "
					  << vd.template getProp<force>(a)[2] <<
	    			  "    DENSITY: " << vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a) <<
					  "    PRESSURE: " << vd.template getProp<Pressure>(a) <<
					  "    POSITION Z: " << vd.getPos(a)[2]
					  << std::endl;*/
	    //	std::cout << "FORCE TERM: " << Tensile(r,rhoa,rhob,Pa,Pb,true) << "      " << massb*Pi(dr,r2,v_rel,rhoa,rhob,massa,massb) << std::endl;
//	    }

	    vd.getPos(a)[0] += dx;
	    vd.getPos(a)[1] += dy;
	    vd.getPos(a)[2] += dz;

	    if (vd.getPos(a)[2] < lower_z)
	    	lower_z = vd.getPos(a)[2];

	    if (vd.getPos(a)[0] <  0.000263878 || vd.getPos(a)[1] < 0.000263878 || vd.getPos(a)[2] < 0.000263878 ||
	    	vd.getPos(a)[0] >  0.000263878+1.59947 || vd.getPos(a)[1] > 0.000263878+0.672972 || vd.getPos(a)[2] > 0.000263878+0.903944)
	    {
//	    	std::cout << "Particle out" << std::endl;
	    	to_remove.add(a.getKey());

	    	particle_out = true;
	    }

	    double velX = vd.template getProp<velocity>(a)[0];
	    double velY = vd.template getProp<velocity>(a)[1];
	    double velZ = vd.template getProp<velocity>(a)[2];
	    double rhop = vd.template getProp<rho>(a);

	    vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<force>(a)[0]*dt2;
	    vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<force>(a)[1]*dt2;
	    vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<force>(a)[2]*dt2;
	    vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);

	    vd.template getProp<velocity_prev>(a)[0] = velX;
	    vd.template getProp<velocity_prev>(a)[1] = velY;
	    vd.template getProp<velocity_prev>(a)[2] = velZ;
	    vd.template getProp<rho_prev>(a) = rhop;

//	    std::cout << "VELOCITY: " << vd.template getProp<velocity>(a)[0] << "     " << vd.template getProp<velocity>(a)[1] << "     " << vd.template getProp<velocity>(a)[2] << std::endl;
//	    std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl;

/*	    if (particle_out == true)
	    {
	    	std::cout << "PARTICLE DENSITY: " << vd.template getProp<rho>(a) << "   PARTICLE PRESSURE: " << vd.template getProp<Pressure>(a) << "     Delta rho: " << vd.template getProp<drho>(a) << "   VELOCITY: " << vd.template getProp<velocity>(a)[0] << "   " << vd.template getProp<velocity>(a)[1] << "    " << vd.template getProp<velocity>(a)[2] << std::endl;
	    }*/

		++part;
	}

	vd.remove(to_remove,0);
}

int main(int argc, char* argv[])
{

	/*!
	 * \page Vector_7_sph_dlb Vector 7 SPH Dam break  simulation with Dynamic load balacing
	 *
	 * ## Initialization ##
	 *
	 * Here we Initialize the library, we create a Box that define our domain, boundary conditions, ghost
	 *
	 * \see \ref e0_s_init
	 *
	 * \snippet Vector/7_SPH_dlb/main.cpp Initialization and parameters
	 *
	 */

	//! \cond [Initialization and parameters] \endcond

    // initialize the library
	openfpm_init(&argc,&argv);

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<3,double> domain({-0.05,-0.05,-0.05},{2.0070,1.0040,1.0040});
	size_t sz[3] = {243,125,125};

	// Fill W_dap
	W_dap = 1.0/Wab(H/1.5);

	// Here we define the boundary conditions of our problem
    size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<3,double> g(2*H);
	
	//! \cond [Initialization and parameters] \endcond

	/*!
	 * \page Vector_7_SPH_dlb Vector 7 SPH Dam break  simulation with Dynamic load balacing
	 *
	 * ## %Vector create ##
	 *
	 * Here we define a distributed vector in 3D, containing 3 properties, a
	 * scalar double, a vector double[3], and a tensor or rank 2 double[3][3].
	 * In this case the vector contain 0 particles initially
	 *
	 * \see \ref vector_inst
	 *
	 * \snippet Vector/1_celllist/main.cpp vector inst
	 *
	 */

	//! \cond [vector inst] \endcond

	particles vd(0,domain,bc,g);

	//! \cond [vector inst] \endcond

	// the scalar is the element at position 0 in the aggregate
	const int type = 0;

//	Box<3,double> fluid_box({dp/2.0,dp/2.0,dp/2.0},{0.4+dp/2.0,0.67-dp/2.0,0.3+dp/2.0});
	Box<3,double> fluid_box({0.2,0.2,0.298},{0.2+2*dp,0.2+2*dp,0.298+2*dp});

	// first we create Fluid particles
	// Fluid particles are created

	auto fluid_it = DrawParticles::DrawBox(vd,sz,domain,fluid_box);
	max_fluid_height = fluid_it.getBoxMargins().getHigh(2);
	h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
	B = (coeff_sound)*(coeff_sound)*gravity*h_swl*rho_zero / gamma_;
	cbar = coeff_sound * sqrt(gravity * h_swl);

//	std::cout << "MAX FLUID: " << max_fluid_height << std::endl;

	while (fluid_it.isNext())
	{
		vd.add();

		vd.getLastPos()[0] = fluid_it.get().get(0);
		vd.getLastPos()[1] = fluid_it.get().get(1);
		vd.getLastPos()[2] = fluid_it.get().get(2);

		vd.template getLastProp<type>() = FLUID;

		// We also initialize the density of the particle and the hydro-static pressure given by
		//
		// rho_zero*g*h = P
		//
		// rho_p = (P/B + 1)^(1/Gamma) * rho_zero
		//

		vd.template getLastProp<Pressure>() = rho_zero * gravity *  (max_fluid_height - fluid_it.get().get(2));

		std::cout << "B: " << B << std::endl;

		vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0/gamma_) * rho_zero;
		vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

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

	// Obstacle

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
	vd.ghost_get<rho,Pressure>();

	auto NN = vd.getCellList(2*H);

	// Evolve


	size_t write = 0;
	size_t it = 0;
	double t = 0.0;
	while (t <= t_end)
	{
		const size_t type = 0;
		const int rho = 1;
		const int Pressure = 2;
		const int drho = 3;
		const int force = 4;
		const int velocity = 5;
		const int velocity_prev = 6;

		vd.map();
		vd.ghost_get<type,rho,Pressure,velocity,velocity_prev>();

		// Calculate pressure from the density
		EqState(vd);

		double max_visc = 0.0;

		// Calc forces
		calc_forces(vd,NN,max_visc);

		// Calculate delta t integration
		double dt = calc_deltaT(vd,max_visc);

//		std::cout << "Calculate deltaT: " << dt << "   " << DtMin << std::endl;

		// VerletStep
		verlet_int(vd,dt);

		t += dt;

		if (write < t*100)
		{

			vd.write("Geometry",write);
			write++;
		}
	}

	//! \cond [finalize] \endcond

	openfpm_finalize();

	//! \cond [finalize] \endcond

	/*!
	 * \page Vector_7_SPH_dlb Vector 7 SPH Dam break  simulation with Dynamic load balacing
	 *
	 * ## Full code ## {#code_e0_sim}
	 *
	 * \include Vector/0_simple/main.cpp
	 *
	 */
}
 
