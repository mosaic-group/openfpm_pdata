#include "Vector/vector_dist.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"
#include "Plot/GoogleChart.hpp"
#include "Plot/util.hpp"
#include "Operators/Vector/vector_dist_operators.hpp"

constexpr int velocity = 0;
constexpr int force = 1;

struct ln_potential
{
	double L;

	ln_potential(double L)
	:L(L){}

	Point<3,double> value(Point<3,double> & xp, Point<3,double> xq, Point<3,double> & sp, Point<3,double> & sq)
	{
		// Energy is a vector we calculate separately Potential and kinetic, and total
		Point<3,double> E;

		float rn = norm(xp - xq) * L;

		E.get(0)= 4.0 * ( 1.0 / (rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn) + 1.0 / ( rn*rn*rn*rn*rn*rn) );
		E.get(1) = (sp * sp) / 2.0;
		E.get(2) = E.get(0) + E.get(1);

		return E;
	}
};

struct ln_force
{
	double L;

	ln_force(double L)
	:L(L){}

	Point<3,double> value(Point<3,double> & xp, Point<3,double> xq)
	{
		Point<3,double> r = xp - xq;

		// take the norm of this vector
		float rn = r.norm();
		r /= rn;
		rn *= L;

		return 24.0*(2.0 / (rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn*rn) -  1.0 / (rn*rn*rn*rn*rn*rn*rn)) * r;
	}
};

void calc_forces(vector_dist<3,double, aggregate<Point<3,double>,Point<3,double>> > & vd, CellList<3, double, FAST, shift<3, double> > & NN, double L)
{

	vd.updateCellList(NN);

	//! \cond [ucl] \endcond

	auto v_force = getV<force>(vd);

	ln_force lf(L);

	v_force = applyKernel_in_sim(v_force,vd,NN,lf);
}

Point<3,double> calc_energy(vector_dist<3,double, aggregate<Point<3,double>,Point<3,double>> > & vd, CellList<3, double, FAST, shift<3, double> > & NN, double L)
{

	Point<3,double> E;

	vd.updateCellList(NN);

	auto v_velocity = getV<velocity>(vd);

	ln_potential lf(L);

	auto eE = applyKernel_reduce(v_velocity,vd,NN,lf);
	eE.init();

	E = eE.value(0);

	return E;
}


int main(int argc, char* argv[])
{
	double dt = 0.005;
	double L = 10.0;
	float r_cut = 0.3;

	openfpm::vector<double> x;
	openfpm::vector<openfpm::vector<double>> y;


	openfpm_init(&argc,&argv);
	Vcluster & v_cl = create_vcluster();

	size_t sz[3] = {10,10,10};
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};
	Ghost<3,float> ghost(r_cut);

	vector_dist<3,double, aggregate<Point<3,double>,Point<3,double>> > vd(0,box,bc,ghost);

	auto v_force = getV<force>(vd);
	auto v_velocity = getV<velocity>(vd);
	auto v_pos = getV<PROP_POS>(vd);

	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		vd.add();

		auto key = it.get();

		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

		++it;
	}

	v_force = 0;
	v_velocity = 0;

	auto NN = vd.getCellList(r_cut);

	// calculate forces
	calc_forces(vd,NN,L);
	unsigned long int f = 0;

	// MD time stepping
	for (size_t i = 0; i < 10000 ; i++)
	{
		auto exp1 = v_velocity + 0.5*dt*v_force;
		auto exp2 = v_pos + v_velocity*dt;

		assign(v_velocity,exp1,
			   v_pos,exp2);

		vd.map();
		vd.template ghost_get<>();

		// calculate forces or a(tn + 1) Step 2
		calc_forces(vd,NN,L);

		v_velocity = v_velocity + 0.5*dt*v_force;

		if (i % 100 == 0)
		{
			vd.deleteGhost();
			vd.write("particles_",f);

			// We calculate the energy
			Point<3,double> energy = calc_energy(vd,NN,L);
			auto & vcl = create_vcluster();
			vcl.sum(energy.get(2));
			vcl.execute();

			// we resync the ghost
			vd.ghost_get<>();

			// we save the energy calculated at time step i c contain the time-step y contain the energy
			x.add(i);
			y.add({energy.get(2)});

			// We also print on terminal the value of the energy
			// only one processor (master) write on terminal
			if (vcl.getProcessUnitID() == 0)
				std::cout << "Energy: " << energy.get(2) << std::endl;

			f++;
		}
	}

	GCoptions options;
	options.title = std::string("Energy with time");
	options.yAxis = std::string("Energy");
	options.xAxis = std::string("iteration");
	options.lineWidth = 1.0;

	GoogleChart cg;
	cg.AddLinesGraph(x,y,options);
	cg.write("gc_plot2_out.html");

	openfpm_finalize();
}
