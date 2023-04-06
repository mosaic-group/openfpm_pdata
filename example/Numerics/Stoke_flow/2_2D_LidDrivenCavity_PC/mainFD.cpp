//
// Created by Abhinav Singh on 15.11.2021.
//

/*!
 * \page Lid_Driven_Cavity_FD Lid driven cavity with Pressure Correction and FD
 *
 * # Lid Driven Cavity Problem with Pressure Correction and FD # {#num_2dlidFD}
 *
 * In this example, we solve the incompressible Navier-Stokes equation in a 2D Square Box:
 *
 * @f[ \mathbf{v}\cdot(\nabla \mathbf{v})-\frac{1}{\text{Re}}\mathrm{\Delta} \mathbf{v}=-\nabla \Pi \label{eq:NS} \\ \nabla\cdot \mathbf{v}=0 \\	\mathbf{v}(x_b,y_b)=(0,0), \text{ except } 	\mathbf{v}(x_b,1)=(1,0)\, , @f]
 *
 * We do that by solving the implicit stokes equation and and employing an iterative pressure correction scheme:
 *
 * Output:
 * Steady State Solution to the Lid Driven Cavity Problem.
 *
 * ## Including the headers ## {#lid_c1_include}
 *
 * These are the header files that we need to include:
 *
 * @snippet Numerics/Stoke_flow/2_2D_LidDrivenCavity_PC/mainFD.cpp LidFDInclude
 *
 */

//! @cond [LidFDInclude] @endcond
// Include Vector Expression,Vector Expressions for Subset,FD , and Solver header files
#include "config.h"
#include <iostream>
#include "FiniteDifference/FD_Solver.hpp"
#include "FiniteDifference/FD_op.hpp"
#include "FiniteDifference/FD_expressions.hpp"
#include "FiniteDifference/util/EqnsStructFD.hpp"
//! \cond [LidFDInclude] \endcond

int main(int argc, char* argv[])
{
    {    /*!
	 * \page Lid_Driven_Cavity_FD Stokes Lid driven cavity with Pressure Correction and FD
	 *
	 * ## Initialization ## {#init2dlidl}
	 *
	 * * Initialize the library
	 * * Define some useful constants
	 * * define Ghost size
	 * * Non-periodic boundary conditions
	 *
     * @snippet Numerics/Stoke_flow/2_2D_LidDrivenCavity_PC/mainFD.cpp LidFDInit
	 *
	 */

	//! \cond [LidFDInit] \endcond
	    openfpm_init(&argc,&argv);
        using namespace FD;
        timer tt2;
        tt2.start();
        size_t gd_sz = 81;
        constexpr int x = 0;
        constexpr int y = 1;
        const size_t szu[2] = {gd_sz,gd_sz};
        int sz[2]={int(gd_sz),int(gd_sz)};
        Box<2, double> box({0, 0}, {1,1});
        periodicity<2> bc = {NON_PERIODIC, NON_PERIODIC};
        double spacing;
        spacing = 1.0 / (sz[0] - 1);
        Ghost<2,long int> ghost(1);
        auto &v_cl = create_vcluster();
        typedef aggregate<double, VectorS<2, double>, VectorS<2, double>,VectorS<2, double>,double,VectorS<2, double>,double,double> LidCavity;
        grid_dist_id<2, double, LidCavity> domain(szu, box,ghost,bc);
        double x0, y0, x1, y1;
        x0 = box.getLow(0);
        y0 = box.getLow(1);
        x1 = box.getHigh(0);
        y1 = box.getHigh(1);
        
      //! \cond [LidFDInit] \endcond

           /*!
	 * \page Lid_Driven_Cavity_FD Stokes Lid driven cavity with Pressure Correction and FD
	 *
	 * ## Creating particles in the 2D domain## {#init2dlidparts}
	 *
     * We set the appropriate subset number 0 for bulk and other for boundary.
     * Note that for different walls we need different subsets as for the pressure, we need normal derivative zero condition.
	 *
     * @snippet Numerics/Stoke_flow/2_2D_LidDrivenCavity_PC/mainFD.cpp LidFDInitPart
	 *
	 */

	//! \cond [LidFDInitPart] \endcond
        auto it = domain.getDomainIterator();
        while (it.isNext())
        {
            auto key = it.get();
            auto gkey = it.getGKey(key);
            double x = gkey.get(0) * domain.spacing(0);
            double y = gkey.get(1) * domain.spacing(1);
            if (y==1)
            {domain.get<1>(key) = 1.0;}
            else
            {domain.get<1>(key) = 0.0;}

            ++it;
        }
        domain.ghost_get<0>();
    //! \cond [LidFDInitPart] \endcond

        /*!
	     * \page Lid_Driven_Cavity_FD Stokes Lid driven cavity with Pressure Correction and FD
         *
         * ##Creating Subsets and Vector Expressions for Fields and Differential Operators for easier encoding. ## {#init2dlidana3}
         *
         * @snippet Numerics/Stoke_flow/2_2D_LidDrivenCavity_PC/mainFD.cpp LidFDexp
         *
         */

	    //! \cond [LidFDexp] \endcond
        Derivative_x Dx;
        Derivative_y Dy;
        Derivative_xx Dxx;
        Derivative_yy Dyy;
        auto P = getV<0>(domain);
        auto V = getV<1>(domain);
        auto RHS = getV<2>(domain);
        auto dV = getV<3>(domain);
        auto div = getV<4>(domain);
        auto V_star = getV<5>(domain);
        auto H = getV<6>(domain);
        auto dP = getV<7>(domain);
        //! \cond [LidFDexp] \endcond

        /*!
	     * \page Lid_Driven_Cavity_FD Stokes Lid driven cavity with Pressure Correction and FD
         *
         * ##Creating a 3D implicit solver for the given set of particles and iteratively solving wit pressure correction. ## {#init3dballana3}
         *
         * @snippet Numerics/Stoke_flow/2_2D_LidDrivenCavity_PC/mainFD.cpp LidFDSol
         *
         */

	    //! \cond [LidFDSol] \endcond
        eq_id x_comp,y_comp;
        x_comp.setId(0);
        y_comp.setId(1);
        double sum, sum1, sum_k,V_err_old;
        auto StokesX=(V[x]*Dx(V_star[x])+V[y]*Dy(V_star[x]))-(1.0/100.0)*(Dxx(V_star[x])+Dyy(V_star[x]));
        auto StokesY=(V[x]*Dx(V_star[y])+V[y]*Dy(V_star[y]))-(1.0/100.0)*(Dxx(V_star[y])+Dyy(V_star[y]));
        RHS=0;
        dV=0;
        double V_err=1;
        int n;
        while (V_err >= 0.05 && n <= 50) {
            if (n%5==0){
                domain.ghost_get<0,1>(SKIP_LABELLING);
                domain.write_frame("LID",n,BINARY);
                domain.ghost_get<0>();
            }
            domain.ghost_get<0>(SKIP_LABELLING);
            RHS[x] = -Dx(P);
            RHS[y] = -Dy(P);
            FD_scheme<equations2d2,decltype(domain)> Solver(ghost,domain);
            Solver.impose(StokesX, {1,1},{sz[0]-2,sz[1]-2}, RHS[x], x_comp);
            Solver.impose(StokesY, {1,1},{sz[0]-2,sz[1]-2}, RHS[y], y_comp);
            Solver.impose(V_star[x], {0,sz[1]-1},{sz[0]-1,sz[1]-1}, 1.0, x_comp);
            Solver.impose(V_star[y], {0,sz[1]-1},{sz[0]-1,sz[1]-1}, 0.0, y_comp);
            Solver.impose(V_star[x], {0,1},{0,sz[1]-2}, 0.0, x_comp);
            Solver.impose(V_star[y], {0,1},{0,sz[1]-2}, 0.0, y_comp);
            Solver.impose(V_star[x], {sz[0]-1,1},{sz[0]-1,sz[1]-2}, 0.0, x_comp);
            Solver.impose(V_star[y], {sz[0]-1,1},{sz[0]-1,sz[1]-2}, 0.0, y_comp);
            Solver.impose(V_star[x], {0,0},{sz[0]-1,0}, 0.0, x_comp);
            Solver.impose(V_star[y], {0,0},{sz[0]-1,0}, 0.0, y_comp);
            Solver.solve(V_star[x], V_star[y]);

            domain.ghost_get<5>(SKIP_LABELLING);
            div = (Dx(V_star[x]) + Dy(V_star[y]));

            FD_scheme<equations2d1E,decltype(domain)> SolverH(ghost,domain);
            auto Helmholtz = Dxx(H)+Dyy(H);
            SolverH.impose(Helmholtz,{1,1},{sz[0]-2,sz[1]-2},div);
            SolverH.impose(Dy(H), {0,sz[0]-1},{sz[0]-1,sz[1]-1},0);
            SolverH.impose(Dx(H), {sz[0]-1,1},{sz[0]-1,sz[1]-2},0);
            SolverH.impose(Dx(H), {0,0},{sz[0]-1,0},0,x_comp,true);
            SolverH.impose(H, {0,0},{0,0},0);
            SolverH.impose(Dy(H), {0,1},{0,sz[1]-2},0);
            SolverH.solve(H);
            domain.ghost_get<1,4,6>(SKIP_LABELLING);
            P = P - 0.01*(div-0.5*(V[x]*Dx(H)+V[y]*Dy(H)));
            V_star[0] = V_star[0] - Dx(H);
            V_star[1] = V_star[1] - Dy(H);
            sum = 0;
            sum1 = 0;
            auto it2 = domain.getDomainIterator();
            while (it2.isNext()) {
                auto p = it2.get();
                sum += (domain.getProp<5>(p)[0] - domain.getProp<1>(p)[0]) *
                       (domain.getProp<5>(p)[0] - domain.getProp<1>(p)[0]) +
                       (domain.getProp<5>(p)[1] - domain.getProp<1>(p)[1]) *
                       (domain.getProp<5>(p)[1] - domain.getProp<1>(p)[1]);
                sum1 += domain.getProp<5>(p)[0] * domain.getProp<5>(p)[0] +
                        domain.getProp<5>(p)[1] * domain.getProp<5>(p)[1];
                ++it2;
            }
            V[x] = V_star[x];
            V[y] = V_star[y];
            v_cl.sum(sum);
            v_cl.sum(sum1);
            v_cl.execute();
            sum = sqrt(sum);
            sum1 = sqrt(sum1);
            V_err_old = V_err;
            V_err = sum / sum1;
            n++;
            if (v_cl.rank() == 0) {
                std::cout << "Rel l2 cgs err in V = " << V_err << "."<< std::endl;
            }
        }
        domain.write("LID_final");
        tt2.stop();
        if (v_cl.rank() == 0) {
            std::cout << "The entire solver  took " << tt2.getcputime() << "(CPU) ------ " << tt2.getwct()
                      << "(Wall) Seconds.";
        }
    }   
    openfpm_finalize();
 //! \cond [LidFDSol] \endcond

         /*!
	     * \page Lid_Driven_Cavity_FD Stokes Lid driven cavity with Pressure Correction and FD
         *
	     * # Full code # {#num_lid_2D_codeFD}
         *
         * \include Numerics/Stoke_flow/2_2D_LidDrivenCavity_PC/mainFD.cpp
         *
         */

}