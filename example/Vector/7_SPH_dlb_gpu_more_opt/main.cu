/*! \page Vector_7_sph_dlb_gpu_more_opt Vector 7 SPH Dam break simulation with Dynamic load balacing on Multi-GPU (more optimized version)
 *
 *
 * [TOC]
 *
 *
 * # SPH with Dynamic load Balancing on GPU (More Optimized) # {#SPH_dlb_gpu_more_opt}
 *
 *
 * This example show the classical SPH Dam break simulation with load balancing and dynamic load balancing. The main difference with
 * \ref{SPH_dlb_gpu_opt} is that here we use 2 kernel to calculate forces one for fluid and one for boundaries. Also we use the function
 * get_indexes_by_type to get the indexes of the fluid and boundary particles and use these two set to launch two distinct kernel
 * (one over fluid and one over boundary) to calculate forces and density change. set. Simulate 1.5 second should be duable on mobile
 *  1050Ti in about 1 hour and 7 minutes
 *
 * \htmlonly
 * <a href="#" onclick="hide_show('vector-video-3')" >Simulation video 1</a><br>
 * <div style="display:none" id="vector-video-3">
 * <video id="vid3" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/7_SPH_dlb/sph_gpu1.mp4" type="video/mp4"></video>
 * </div>
 * <a href="#" onclick="hide_show('vector-video-4')" >Simulation video 2</a><br>
 * <div style="display:none" id="vector-video-4">
 * <video id="vid4" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/7_SPH_dlb/sph_gpu2.mp4" type="video/mp4"></video>
 * </div>
 * <a href="#" onclick="hide_show('vector-video-15')" >Simulation video 3</a><br>
 * <div style="display:none" id="vector-video-15">
 * <video id="vid15" width="1200" height="576" controls> <source src="http://openfpm.mpi-cbg.de/web/images/examples/7_SPH_dlb/sph_gpu3.mp4" type="video/mp4"></video>
 * </div>
 * \endhtmlonly
 *
 *
 * ## get_indexes_by_type ## {#e7_sph_more_opt_gibt}
 *
 * This function can be used to get the indexes of a certain type on a particle set and save such indexes in an openfpm::vector<aggregate<unsigned int>>
 * the constructed set of indices can be used to run a kernel on a specific set of particles.
 *
 * \snippet Vector/7_SPH_dlb_gpu_more_opt/main.cu get indexes by type
 *
 * the function get_indexes_by_type has three arguments the first is the vector of the properties of the particles. In
 * this case because we use the sorted particles to calculate forces, so we have to get the indexes for the sorted
 * particles with vd.getPropVectorSort(). In case we want to use the non sorted we use vd.getPropVector(). The second
 * argument is the output containing the indexes of the particles types we want to get. Because the vector can contain
 * ghost particles and real particles setting with the third argument we indicate we want only real particles and no ghost particles
 * The last argument is the GPU context handle
 *
 * we report the full code here
 *
 *
 */

#ifdef __NVCC__

#define PRINT_STACKTRACE
#define STOP_ON_ERROR
#define OPENMPI
#define SCAN_WITH_CUB
#define SORT_WITH_CUB
//#define SE_CLASS1

//#define USE_LOW_REGISTER_ITERATOR

#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"



typedef float real_number;

// A constant to indicate boundary particles
#define BOUNDARY 0

// A constant to indicate fluid particles
#define FLUID 1

// initial spacing between particles dp in the formulas
const real_number dp = 0.00425;
// Maximum height of the fluid water
// is going to be calculated and filled later on
real_number h_swl = 0.0;

// c_s in the formulas (constant used to calculate the sound speed)
const real_number coeff_sound = 20.0;

// gamma in the formulas
const real_number gamma_ = 7.0;

// sqrt(3.0*dp*dp) support of the kernel
const real_number H = 0.00736121593217;

// Eta in the formulas
const real_number Eta2 = 0.01 * H*H;

const real_number FourH2 = 4.0 * H*H;

// alpha in the formula
const real_number visco = 0.1;

// cbar in the formula (calculated later)
real_number cbar = 0.0;

// Mass of the fluid particles
const real_number MassFluid = 0.0000767656;

// Mass of the boundary particles
const real_number MassBound = 0.0000767656;

//

// End simulation time
#ifdef TEST_RUN
const real_number t_end = 0.001;
#else
const real_number t_end = 1.5;
#endif

// Gravity acceleration
const real_number gravity = 9.81;

// Reference densitu 1000Kg/m^3
const real_number rho_zero = 1000.0;

// Filled later require h_swl, it is b in the formulas
real_number B = 0.0;

// Constant used to define time integration
const real_number CFLnumber = 0.2;

// Minimum T
const real_number DtMin = 0.00001;

// Minimum Rho allowed
const real_number RhoMin = 700.0;

// Maximum Rho allowed
const real_number RhoMax = 1300.0;

// Filled in initialization
real_number max_fluid_height = 0.0;

// Properties

// FLUID or BOUNDARY
const size_t type = 0;

// Density
const int rho = 1;

// Density at step n-1
const int rho_prev = 2;

// Pressure
const int Pressure = 3;

// Delta rho calculated in the force calculation
const int drho = 4;

// calculated force
const int force = 5;

// velocity
const int velocity = 6;

// velocity at previous step
const int velocity_prev = 7;

const int red = 8;

const int red2 = 9;

// Type of the vector containing particles
typedef vector_dist_gpu<3,real_number,aggregate<unsigned int,real_number,  real_number,    real_number,     real_number,     real_number[3], real_number[3], real_number[3], real_number, real_number>> particles;
//                                              |          |             |               |                |                |               |               |               |            |
//                                              |          |             |               |                |                |               |               |               |            |
//                                             type      density       density        Pressure          delta            force          velocity        velocity        reduction     another
//                                                                     at n-1                           density                                         at n - 1        buffer        reduction buffer


struct ModelCustom
{
	template<typename Decomposition, typename vector> inline void addComputation(Decomposition & dec,
			                                                                     vector & vd,
																				 size_t v,
																				 size_t p)
	{
		if (vd.template getProp<type>(p) == FLUID)
			dec.addComputationCost(v,4);
		else
			dec.addComputationCost(v,3);
	}

	template<typename Decomposition> inline void applyModel(Decomposition & dec, size_t v)
	{
		dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
	}

	real_number distributionTol()
	{
		return 1.01;
	}
};

template<typename vd_type>
__global__ void EqState_gpu(vd_type vd, real_number B)
{
	auto a = GET_PARTICLE(vd);

	real_number rho_a = vd.template getProp<rho>(a);
	real_number rho_frac = rho_a / rho_zero;

	vd.template getProp<Pressure>(a) = B*( rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac*rho_frac - 1.0);
}

inline void EqState(particles & vd)
{
	auto it = vd.getDomainIteratorGPU();

	CUDA_LAUNCH(EqState_gpu,it,vd.toKernel(),B);
}


const real_number a2 = 1.0/M_PI/H/H/H;

inline __device__ __host__ real_number Wab(real_number r)
{
	r /= H;

	if (r < 1.0)
		return (1.0 - 3.0/2.0*r*r + 3.0/4.0*r*r*r)*a2;
	else if (r < 2.0)
		return (1.0/4.0*(2.0 - r*r)*(2.0 - r*r)*(2.0 - r*r))*a2;
	else
		return 0.0;
}


const real_number c1 = -3.0/M_PI/H/H/H/H;
const real_number d1 = 9.0/4.0/M_PI/H/H/H/H;
const real_number c2 = -3.0/4.0/M_PI/H/H/H/H;
const real_number a2_4 = 0.25*a2;
// Filled later
real_number W_dap = 0.0;

inline __device__ __host__ void DWab(Point<3,real_number> & dx, Point<3,real_number> & DW, real_number r)
{
	const real_number qq=r/H;

    real_number qq2 = qq * qq;
    real_number fac1 = (c1*qq + d1*qq2)/r;
    real_number b1 = (qq < 1.0f)?1.0f:0.0f;

    real_number wqq = (2.0f - qq);
    real_number fac2 = c2 * wqq * wqq / r;
    real_number b2 = (qq >= 1.0f && qq < 2.0f)?1.0f:0.0f;

    real_number factor = (b1*fac1 + b2*fac2);

    DW.get(0) = factor * dx.get(0);
    DW.get(1) = factor * dx.get(1);
    DW.get(2) = factor * dx.get(2);
}

// Tensile correction
inline __device__ __host__  real_number Tensile(real_number r, real_number rhoa, real_number rhob, real_number prs1, real_number prs2, real_number W_dap)
{
	const real_number qq=r/H;
	//-Cubic Spline kernel
	real_number wab;
	if(r>H)
	{
		real_number wqq1=2.0f-qq;
		real_number wqq2=wqq1*wqq1;

		wab=a2_4*(wqq2*wqq1);
	}
	else
	{
	    real_number wqq2=qq*qq;
	    real_number wqq3=wqq2*qq;

	    wab=a2*(1.0f-1.5f*wqq2+0.75f*wqq3);
	}

	//-Tensile correction.
	real_number fab=wab*W_dap;
	fab*=fab; fab*=fab; //fab=fab^4
	const real_number tensilp1=(prs1/(rhoa*rhoa))*(prs1>0.0f? 0.01f: -0.2f);
	const real_number tensilp2=(prs2/(rhob*rhob))*(prs2>0.0f? 0.01f: -0.2f);

	return (fab*(tensilp1+tensilp2));
}


inline __device__ __host__ real_number Pi(const Point<3,real_number> & dr, real_number rr2, Point<3,real_number> & dv, real_number rhoa, real_number rhob, real_number massb, real_number cbar, real_number & visc)
{
	const real_number dot = dr.get(0)*dv.get(0) + dr.get(1)*dv.get(1) + dr.get(2)*dv.get(2);
	const real_number dot_rr2 = dot/(rr2+Eta2);
	visc=(dot_rr2 < visc)?visc:dot_rr2;

	if(dot < 0)
	{
		const float amubar=H*dot_rr2;
		const float robar=(rhoa+rhob)*0.5f;
		const float pi_visc=(-visco*cbar*amubar/robar);

		return pi_visc;
    }
	else
		return 0.0f;
}

template<typename particles_type, typename fluid_ids_type,typename NN_type>
__global__ void calc_forces_fluid_gpu(particles_type vd, fluid_ids_type fids, NN_type NN, real_number W_dap, real_number cbar)
{
	// ... a
	unsigned int a;

	GET_PARTICLE_BY_ID(a,fids);

	real_number max_visc = 0.0f;

	// Get the position xp of the particle
	Point<3,real_number> xa = vd.getPos(a);

	// Type of the particle
	unsigned int typea = vd.template getProp<type>(a);

	// Get the density of the of the particle a
	real_number rhoa = vd.template getProp<rho>(a);

	// Get the pressure of the particle a
	real_number Pa = vd.template getProp<Pressure>(a);

	// Get the Velocity of the particle a
	Point<3,real_number> va = vd.template getProp<velocity>(a);

	Point<3,real_number> force_;
	force_.get(0) = 0.0f;
	force_.get(1) = 0.0f;
	force_.get(2) = -gravity;
	real_number drho_ = 0.0f;

	// Get an iterator over the neighborhood particles of p
	auto Np = NN.getNNIteratorBox(NN.getCell(xa));

	// For each neighborhood particle
	while (Np.isNext() == true)
	{
		// ... q
		auto b = Np.get_sort();

		// Get the position xp of the particle
		Point<3,real_number> xb = vd.getPos(b);

		// if (p == q) skip this particle this condition should be done in the r^2 = 0
		//if (a == b)	{++Np; continue;};

        unsigned int typeb = vd.template getProp<type>(b);

        real_number massb = (typeb == FLUID)?MassFluid:MassBound;
        Point<3,real_number> vb = vd.template getProp<velocity>(b);
        real_number Pb = vd.template getProp<Pressure>(b);
        real_number rhob = vd.template getProp<rho>(b);

		// Get the distance between p and q
		Point<3,real_number> dr = xa - xb;
		Point<3,real_number> v_rel = va - vb;
		// take the norm of this vector
		real_number r2 = norm2(dr);

		// if they interact
		if (r2 < FourH2 && r2 >= 1e-16)
		{
			real_number r = sqrtf(r2);

			Point<3,real_number> DW;
			DWab(dr,DW,r);

			real_number factor = - massb*((Pa + Pb) / (rhoa * rhob) + Tensile(r,rhoa,rhob,Pa,Pb,W_dap) + Pi(dr,r2,v_rel,rhoa,rhob,massb,cbar,max_visc));

			// Bound - Bound does not produce any change
			factor = (typea == BOUNDARY && typeb == BOUNDARY)?0.0f:factor;

			force_.get(0) += factor * DW.get(0);
			force_.get(1) += factor * DW.get(1);
			force_.get(2) += factor * DW.get(2);

			real_number scal = massb*(v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2));
			scal = (typea == BOUNDARY && typeb == BOUNDARY)?0.0f:scal;

			drho_ += scal;
		}

		++Np;
	}

	vd.template getProp<red>(a) = max_visc;

	vd.template getProp<force>(a)[0] = force_.get(0);
	vd.template getProp<force>(a)[1] = force_.get(1);
	vd.template getProp<force>(a)[2] = force_.get(2);
	vd.template getProp<drho>(a) = drho_;
}

template<typename particles_type, typename fluid_ids_type,typename NN_type>
__global__ void calc_forces_border_gpu(particles_type vd, fluid_ids_type fbord, NN_type NN, real_number W_dap, real_number cbar)
{
	// ... a
	unsigned int a;

	GET_PARTICLE_BY_ID(a,fbord);

	real_number max_visc = 0.0f;

	// Get the position xp of the particle
	Point<3,real_number> xa = vd.getPos(a);

	// Type of the particle
	unsigned int typea = vd.template getProp<type>(a);

	// Get the Velocity of the particle a
	Point<3,real_number> va = vd.template getProp<velocity>(a);

	real_number drho_ = 0.0f;

	// Get an iterator over the neighborhood particles of p
	auto Np = NN.getNNIteratorBox(NN.getCell(xa));

	// For each neighborhood particle
	while (Np.isNext() == true)
	{
		// ... q
		auto b = Np.get_sort();

		// Get the position xp of the particle
		Point<3,real_number> xb = vd.getPos(b);

		// if (p == q) skip this particle this condition should be done in the r^2 = 0
		//if (a == b)	{++Np; continue;};

        unsigned int typeb = vd.template getProp<type>(b);

        real_number massb = (typeb == FLUID)?MassFluid:MassBound;
        Point<3,real_number> vb = vd.template getProp<velocity>(b);

		// Get the distance between p and q
		Point<3,real_number> dr = xa - xb;
		Point<3,real_number> v_rel = va - vb;
		// take the norm of this vector
		real_number r2 = norm2(dr);

		// if they interact
		if (r2 < FourH2 && r2 >= 1e-16)
		{
			real_number r = sqrtf(r2);

			Point<3,real_number> DW;
			DWab(dr,DW,r);

			real_number scal = massb*(v_rel.get(0)*DW.get(0)+v_rel.get(1)*DW.get(1)+v_rel.get(2)*DW.get(2));
			scal = (typea == BOUNDARY && typeb == BOUNDARY)?0.0f:scal;

			drho_ += scal;
		}

		++Np;
	}

	vd.template getProp<red>(a) = max_visc;

	vd.template getProp<drho>(a) = drho_;
}

struct type_is_fluid
{
	__device__ static bool check(int c)
	{
		return c == FLUID;
	}
};

struct type_is_border
{
        __device__ static bool check(int c)
        {
                return c == BOUNDARY;
        }
};

template<typename CellList> inline void calc_forces(particles & vd, CellList & NN, real_number & max_visc, size_t cnt, openfpm::vector_gpu<aggregate<int>> & fluid_ids, openfpm::vector_gpu<aggregate<int>> & border_ids)
{
	// Update the cell-list
	vd.updateCellList<type,rho,Pressure,velocity>(NN);

	//! \cond [get indexes by type] \endcond

	// get the particles fluid ids
	get_indexes_by_type<type,type_is_fluid>(vd.getPropVectorSort(),fluid_ids,vd.size_local(),vd.getVC().getGpuContext());

	// get the particles fluid ids
	get_indexes_by_type<type,type_is_border>(vd.getPropVectorSort(),border_ids,vd.size_local(),vd.getVC().getGpuContext());

	auto part = fluid_ids.getGPUIterator(96);
	CUDA_LAUNCH(calc_forces_fluid_gpu,part,vd.toKernel_sorted(),fluid_ids.toKernel(),NN.toKernel(),W_dap,cbar);

	part = border_ids.getGPUIterator(96);
	CUDA_LAUNCH(calc_forces_border_gpu,part,vd.toKernel_sorted(),border_ids.toKernel(),NN.toKernel(),W_dap,cbar);

	//! \cond [get indexes by type] \endcond

	vd.merge_sort<force,drho,red>(NN);

	max_visc = reduce_local<red,_max_>(vd);
}

template<typename vector_type>
__global__ void max_acceleration_and_velocity_gpu(vector_type vd)
{
	auto a = GET_PARTICLE(vd);

	Point<3,real_number> acc(vd.template getProp<force>(a));
	vd.template getProp<red>(a) = norm(acc);

	Point<3,real_number> vel(vd.template getProp<velocity>(a));
	vd.template getProp<red2>(a) = norm(vel);
}

void max_acceleration_and_velocity(particles & vd, real_number & max_acc, real_number & max_vel)
{
	// Calculate the maximum acceleration
	auto part = vd.getDomainIteratorGPU();

	CUDA_LAUNCH(max_acceleration_and_velocity_gpu,part,vd.toKernel());

	max_acc = reduce_local<red,_max_>(vd);
	max_vel = reduce_local<red2,_max_>(vd);

	Vcluster<> & v_cl = create_vcluster();
	v_cl.max(max_acc);
	v_cl.max(max_vel);
	v_cl.execute();
}


real_number calc_deltaT(particles & vd, real_number ViscDtMax)
{
	real_number Maxacc = 0.0;
	real_number Maxvel = 0.0;
	max_acceleration_and_velocity(vd,Maxacc,Maxvel);

	//-dt1 depends on force per unit mass.
	const real_number dt_f = (Maxacc)?sqrt(H/Maxacc):std::numeric_limits<float>::max();

	//-dt2 combines the Courant and the viscous time-step controls.
	const real_number dt_cv = H/(std::max(cbar,Maxvel*10.f) + H*ViscDtMax);

	//-dt new value of time step.
	real_number dt=real_number(CFLnumber)*std::min(dt_f,dt_cv);
	if(dt<real_number(DtMin))
	{dt=real_number(DtMin);}

	return dt;
}

template<typename vector_dist_type>
__global__ void verlet_int_gpu(vector_dist_type vd, real_number dt, real_number dt2, real_number dt205)
{
	// ... a
	auto a = GET_PARTICLE(vd);

	// if the particle is boundary
	if (vd.template getProp<type>(a) == BOUNDARY)
	{
		// Update rho
		real_number rhop = vd.template getProp<rho>(a);

		// Update only the density
    	vd.template getProp<velocity>(a)[0] = 0.0;
    	vd.template getProp<velocity>(a)[1] = 0.0;
    	vd.template getProp<velocity>(a)[2] = 0.0;
    	real_number rhonew = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);
    	vd.template getProp<rho>(a) = (rhonew < rho_zero)?rho_zero:rhonew;

	    vd.template getProp<rho_prev>(a) = rhop;

	    vd.template getProp<red>(a) = 0;

		return;
	}

	//-Calculate displacement and update position
	real_number dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
    real_number dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
    real_number dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;

    vd.getPos(a)[0] += dx;
    vd.getPos(a)[1] += dy;
    vd.getPos(a)[2] += dz;

    real_number velX = vd.template getProp<velocity>(a)[0];
    real_number velY = vd.template getProp<velocity>(a)[1];
    real_number velZ = vd.template getProp<velocity>(a)[2];

    real_number rhop = vd.template getProp<rho>(a);

	vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<force>(a)[0]*dt2;
	vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<force>(a)[1]*dt2;
	vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<force>(a)[2]*dt2;
	vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2*vd.template getProp<drho>(a);

    // Check if the particle go out of range in space and in density
    if (vd.getPos(a)[0] <  0.0 || vd.getPos(a)[1] < 0.0 || vd.getPos(a)[2] < 0.0 ||
        vd.getPos(a)[0] >  1.61 || vd.getPos(a)[1] > 0.68 || vd.getPos(a)[2] > 0.50 ||
		vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
    {
    	vd.template getProp<red>(a) = 1;
    }
    else
    {
    	vd.template getProp<red>(a) = 0;
    }


    vd.template getProp<velocity_prev>(a)[0] = velX;
    vd.template getProp<velocity_prev>(a)[1] = velY;
    vd.template getProp<velocity_prev>(a)[2] = velZ;
    vd.template getProp<rho_prev>(a) = rhop;
}

size_t cnt = 0;

void verlet_int(particles & vd, real_number dt)
{
	// particle iterator
	auto part = vd.getDomainIteratorGPU();

	real_number dt205 = dt*dt*0.5;
	real_number dt2 = dt*2.0;

	CUDA_LAUNCH(verlet_int_gpu,part,vd.toKernel(),dt,dt2,dt205);

	// remove the particles marked
	remove_marked<red>(vd);

	// increment the iteration counter
	cnt++;
}

template<typename vector_type>
__global__ void euler_int_gpu(vector_type vd,real_number dt, real_number dt205)
{
	// ... a
	auto a = GET_PARTICLE(vd);

	// if the particle is boundary
	if (vd.template getProp<type>(a) == BOUNDARY)
	{
		// Update rho
		real_number rhop = vd.template getProp<rho>(a);

		// Update only the density
    	vd.template getProp<velocity>(a)[0] = 0.0;
    	vd.template getProp<velocity>(a)[1] = 0.0;
    	vd.template getProp<velocity>(a)[2] = 0.0;
    	real_number rhonew = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);
    	vd.template getProp<rho>(a) = (rhonew < rho_zero)?rho_zero:rhonew;

	    vd.template getProp<rho_prev>(a) = rhop;

	    vd.template getProp<red>(a) = 0;

		return;
	}

	//-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
	real_number dx = vd.template getProp<velocity>(a)[0]*dt + vd.template getProp<force>(a)[0]*dt205;
    real_number dy = vd.template getProp<velocity>(a)[1]*dt + vd.template getProp<force>(a)[1]*dt205;
    real_number dz = vd.template getProp<velocity>(a)[2]*dt + vd.template getProp<force>(a)[2]*dt205;

    vd.getPos(a)[0] += dx;
    vd.getPos(a)[1] += dy;
    vd.getPos(a)[2] += dz;

    real_number velX = vd.template getProp<velocity>(a)[0];
    real_number velY = vd.template getProp<velocity>(a)[1];
    real_number velZ = vd.template getProp<velocity>(a)[2];
    real_number rhop = vd.template getProp<rho>(a);

	vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + vd.template getProp<force>(a)[0]*dt;
	vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + vd.template getProp<force>(a)[1]*dt;
   	vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + vd.template getProp<force>(a)[2]*dt;
   	vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt*vd.template getProp<drho>(a);

    // Check if the particle go out of range in space and in density
    if (vd.getPos(a)[0] <  0.0 || vd.getPos(a)[1] < 0.0 || vd.getPos(a)[2] < 0.0 ||
        vd.getPos(a)[0] >  1.61 || vd.getPos(a)[1] > 0.68 || vd.getPos(a)[2] > 0.50 ||
		vd.template getProp<rho>(a) < RhoMin || vd.template getProp<rho>(a) > RhoMax)
    {vd.template getProp<red>(a) = 1;}
    else
    {vd.template getProp<red>(a) = 0;}

    vd.template getProp<velocity_prev>(a)[0] = velX;
    vd.template getProp<velocity_prev>(a)[1] = velY;
    vd.template getProp<velocity_prev>(a)[2] = velZ;
    vd.template getProp<rho_prev>(a) = rhop;
}

void euler_int(particles & vd, real_number dt)
{

	// particle iterator
	auto part = vd.getDomainIteratorGPU();

	real_number dt205 = dt*dt*0.5;

	CUDA_LAUNCH(euler_int_gpu,part,vd.toKernel(),dt,dt205);

	// remove the particles
	remove_marked<red>(vd);

	cnt++;
}

template<typename vector_type, typename NN_type>
__global__ void sensor_pressure_gpu(vector_type vd, NN_type NN, Point<3,real_number> probe, real_number * press_tmp)
{
	real_number tot_ker = 0.0;

	// Get the position of the probe i
	Point<3,real_number> xp = probe;

	// get the iterator over the neighbohood particles of the probes position
	auto itg = NN.getNNIteratorBox(NN.getCell(xp));
	while (itg.isNext())
	{
		auto q = itg.get_sort();

		// Only the fluid particles are importants
		if (vd.template getProp<type>(q) != FLUID)
		{
			++itg;
			continue;
		}

		// Get the position of the neighborhood particle q
		Point<3,real_number> xq = vd.getPos(q);

		// Calculate the contribution of the particle to the pressure
		// of the probe
		real_number r = sqrt(norm2(xp - xq));

		real_number ker = Wab(r) * (MassFluid / rho_zero);

		// Also keep track of the calculation of the summed
		// kernel
		tot_ker += ker;

		// Add the total pressure contribution
		*press_tmp += vd.template getProp<Pressure>(q) * ker;

		// next neighborhood particle
		++itg;
	}

	// We calculate the pressure normalizing the
	// sum over all kernels
	if (tot_ker == 0.0)
	{*press_tmp = 0.0;}
	else
	{*press_tmp = 1.0 / tot_ker * *press_tmp;}
}

template<typename Vector, typename CellList>
inline void sensor_pressure(Vector & vd,
                            CellList & NN,
                            openfpm::vector<openfpm::vector<real_number>> & press_t,
                            openfpm::vector<Point<3,real_number>> & probes)
{
    Vcluster<> & v_cl = create_vcluster();

    press_t.add();

    for (size_t i = 0 ; i < probes.size() ; i++)
    {
    	// A float variable to calculate the pressure of the problem
    	CudaMemory press_tmp_(sizeof(real_number));
    	real_number press_tmp;

        // if the probe is inside the processor domain
		if (vd.getDecomposition().isLocal(probes.get(i)) == true)
		{
			CUDA_LAUNCH_DIM3(sensor_pressure_gpu,1,1,vd.toKernel_sorted(),NN.toKernel(),probes.get(i),(real_number *)press_tmp_.toKernel());

			vd.merge<Pressure>(NN);

			// move calculated pressure on
			press_tmp_.deviceToHost();
			press_tmp = *(real_number *)press_tmp_.getPointer();
		}

		// This is not necessary in principle, but if you
		// want to make all processor aware of the history of the calculated
		// pressure we have to execute this
		v_cl.sum(press_tmp);
		v_cl.execute();

		// We add the calculated pressure into the history
		press_t.last().add(press_tmp);
	}
}

int main(int argc, char* argv[])
{
    // initialize the library
	openfpm_init(&argc,&argv);

	openfpm::vector_gpu<aggregate<int>> fluid_ids;
	openfpm::vector_gpu<aggregate<int>> border_ids;

#ifdef CUDIFY_USE_CUDA
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
#endif

	// It contain for each time-step the value detected by the probes
	openfpm::vector<openfpm::vector<real_number>> press_t;
	openfpm::vector<Point<3,real_number>> probes;

	probes.add({0.8779,0.3,0.02});
	probes.add({0.754,0.31,0.02});

	// Here we define our domain a 2D box with internals from 0 to 1.0 for x and y
	Box<3,real_number> domain({-0.05,-0.05,-0.05},{1.7010,0.7065,0.511});
	size_t sz[3] = {413,179,133};

	// Fill W_dap
	W_dap = 1.0/Wab(H/1.5);

	// Here we define the boundary conditions of our problem
    size_t bc[3]={NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// extended boundary around the domain, and the processor domain
	Ghost<3,real_number> g(2*H);

	particles vd(0,domain,bc,g,DEC_GRAN(128));

	//! \cond [draw fluid] \endcond

	// You can ignore all these dp/2.0 is a trick to reach the same initialization
	// of Dual-SPH that use a different criteria to draw particles
	Box<3,real_number> fluid_box({dp/2.0,dp/2.0,dp/2.0},{0.4+dp/2.0,0.67-dp/2.0,0.3+dp/2.0});

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
	Box<3,real_number> recipient1({0.0,0.0,0.0},{1.6+dp/2.0,0.67+dp/2.0,0.4+dp/2.0});
	Box<3,real_number> recipient2({dp,dp,dp},{1.6-dp/2.0,0.67-dp/2.0,0.4+dp/2.0});

	Box<3,real_number> obstacle1({0.9,0.24-dp/2.0,0.0},{1.02+dp/2.0,0.36,0.45+dp/2.0});
	Box<3,real_number> obstacle2({0.9+dp,0.24+dp/2.0,0.0},{1.02-dp/2.0,0.36-dp,0.45-dp/2.0});
	Box<3,real_number> obstacle3({0.9+dp,0.24,0.0},{1.02,0.36,0.45});

	openfpm::vector<Box<3,real_number>> holes;
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

	// Now that we fill the vector with particles
	ModelCustom md;

	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map();

	///////////////////////////

	// Ok the initialization is done on CPU on GPU we are doing the main loop, so first we offload all properties on GPU

	vd.hostToDevicePos();
	vd.template hostToDeviceProp<type,rho,rho_prev,Pressure,velocity>();


	vd.ghost_get<type,rho,Pressure,velocity>(RUN_ON_DEVICE);

	auto NN = vd.getCellListGPU/*<CELLLIST_GPU_SPARSE<3,float>>*/(2*H / 2.0);
	NN.setBoxNN(2);

	timer tot_sim;
	tot_sim.start();

	size_t write = 0;
	size_t it = 0;
	size_t it_reb = 0;
	real_number t = 0.0;
	while (t <= t_end)
	{
		Vcluster<> & v_cl = create_vcluster();
		timer it_time;
		it_time.start();

		////// Do rebalancing every 200 timesteps
		it_reb++;
		if (it_reb == 300)
		{
			vd.map(RUN_ON_DEVICE);

			// Rebalancer for now work on CPU , so move to CPU
            vd.deviceToHostPos();
            vd.template deviceToHostProp<type>();

			it_reb = 0;
			ModelCustom md;
			vd.addComputationCosts(md);
			vd.getDecomposition().decompose();

			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "REBALANCED " << it_reb << std::endl;}
		}

		vd.map(RUN_ON_DEVICE);

		// Calculate pressure from the density
		EqState(vd);

		real_number max_visc = 0.0;

		vd.ghost_get<type,rho,Pressure,velocity>(RUN_ON_DEVICE);


		// Calc forces
		calc_forces(vd,NN,max_visc,cnt,fluid_ids,border_ids);

		// Get the maximum viscosity term across processors
		v_cl.max(max_visc);
		v_cl.execute();

		// Calculate delta t integration
		real_number dt = calc_deltaT(vd,max_visc);

		// VerletStep or euler step
		it++;
		if (it < 40)
			verlet_int(vd,dt);
		else
		{
			euler_int(vd,dt);
			it = 0;
		}

		t += dt;

		if (write < t*10)
		{
			// Sensor pressure require update ghost, so we ensure that particles are distributed correctly
			// and ghost are updated
			vd.map(RUN_ON_DEVICE);
			vd.ghost_get<type,rho,Pressure,velocity>(RUN_ON_DEVICE);
			vd.updateCellList(NN);

			// calculate the pressure at the sensor points
			//sensor_pressure(vd,NN,press_t,probes);

			std::cout << "OUTPUT " << dt << std::endl;

			// When we write we have move all the particles information back to CPU

			vd.deviceToHostPos();
			vd.deviceToHostProp<type,rho,rho_prev,Pressure,drho,force,velocity,velocity_prev,red,red2>();

			// We copy on another vector with less properties to reduce the size of the output
			vector_dist_gpu<3,real_number,aggregate<unsigned int,real_number[3]>> vd_out(vd.getDecomposition(),0);

			auto ito = vd.getDomainIterator();

			while(ito.isNext())
			{
				auto p = ito.get();

				vd_out.add();

				vd_out.getLastPos()[0] = vd.getPos(p)[0];
				vd_out.getLastPos()[1] = vd.getPos(p)[1];
				vd_out.getLastPos()[2] = vd.getPos(p)[2];

				vd_out.template getLastProp<0>() = vd.template getProp<type>(p);

				vd_out.template getLastProp<1>()[0] = vd.template getProp<velocity>(p)[0];
				vd_out.template getLastProp<1>()[1] = vd.template getProp<velocity>(p)[1];
				vd_out.template getLastProp<1>()[2] = vd.template getProp<velocity>(p)[2];

				++ito;
			}

			vd_out.write_frame("Particles",write,VTK_WRITER | FORMAT_BINARY);
			write++;

			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << it_reb << "   " << cnt << " Max visc: " << max_visc << "   " << vd.size_local()  << std::endl;}
		}
		else
		{
			if (v_cl.getProcessUnitID() == 0)
			{std::cout << "TIME: " << t << "  " << it_time.getwct() << "   " << it_reb << "   " << cnt  << " Max visc: " << max_visc << "   " << vd.size_local() << std::endl;}
		}
	}

	tot_sim.stop();
	std::cout << "Time to complete: " << tot_sim.getwct() << " seconds" << std::endl;


	openfpm_finalize();
}
 
#else

int main(int argc, char* argv[])
{
        return 0;
}

#endif
