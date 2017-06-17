/*
 * CG.hpp
 *
 *  Created on: Jun 16, 2017
 *      Author: i-bird
 */

#ifndef EXAMPLE_NUMERICS_VORTEX_IN_CELL_CG_HPP_
#define EXAMPLE_NUMERICS_VORTEX_IN_CELL_CG_HPP_

typedef grid_dist_id<3,float,aggregate<float>> grid_type_scal;

float norm(grid_type_scal & scl)
{
	double nrm = 0.0;

	Vcluster & v_cl = create_vcluster();

	auto it = scl.getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();

		nrm += scl.template getProp<0>(p)*scl.template getProp<0>(p);

		++it;
	}

	v_cl.sum(nrm);
	v_cl.execute();

	return nrm;
}

void copy(grid_type_scal & src, grid_type_scal & dst)
{
	auto it = src.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		dst.getProp<0>(key) = src.getProp<0>(key);

		++it;
	}
}


void x_plus_alpha_p(grid_type_scal & x, grid_type_scal & p, float alpha)
{
	auto it = x.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		x.getProp<0>(key) = x.getProp<0>(key) + alpha * p.getProp<0>(key);

		++it;
	}
}

void x_plus_alpha_p(grid_type_scal & xn, grid_type_scal & x, grid_type_scal & p, float alpha)
{
	auto it = x.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		xn.getProp<0>(key) = x.getProp<0>(key) + alpha * p.getProp<0>(key);

		++it;
	}
}

float alpha(void (* A)(grid_type_scal & in, grid_type_scal & out),
		   grid_type_scal & r,
		   grid_type_scal & p,
		   grid_type_scal & Ap)
{
	auto it = r.getDomainIterator();
	double scal_r = 0.0;
	double scal_pAp = 0.0;

	A(p,Ap);

	while (it.isNext())
	{
		auto key = it.get();

		scal_r += r.template getProp<0>(key)*r.template getProp<0>(key);
		scal_pAp += p.template getProp<0>(key)*Ap.template getProp<0>(key);


		++it;
	}

	Vcluster & v_cl = create_vcluster();
	v_cl.sum(scal_r);
	v_cl.sum(scal_pAp);
	v_cl.execute();

	return scal_r / scal_pAp;
}

float beta(grid_type_scal & r,
		   grid_type_scal & rn)
{
	auto it = r.getDomainIterator();
	double scal_r = 0.0;
	double scal_rn = 0.0;

	while (it.isNext())
	{
		auto key = it.get();

		scal_r += r.template getProp<0>(key)*r.template getProp<0>(key);
		scal_rn += rn.template getProp<0>(key)*rn.template getProp<0>(key);

		++it;
	}

	Vcluster & v_cl = create_vcluster();
	v_cl.sum(scal_r);
	v_cl.sum(scal_rn);
	v_cl.execute();

	return scal_rn / scal_r;
}

void CG(void (* A)(grid_type_scal & in, grid_type_scal & out), grid_type_scal & x, grid_type_scal & b)
{
	grid_type_scal tmp(b.getDecomposition(),b.getGridInfo().getSize(),b.getDecomposition().getGhost());
	grid_type_scal r(b.getDecomposition(),b.getGridInfo().getSize(),b.getDecomposition().getGhost());
	grid_type_scal rn(b.getDecomposition(),b.getGridInfo().getSize(),b.getDecomposition().getGhost());
	grid_type_scal p(b.getDecomposition(),b.getGridInfo().getSize(),b.getDecomposition().getGhost());
	grid_type_scal Ap(b.getDecomposition(),b.getGridInfo().getSize(),b.getDecomposition().getGhost());


	// r_0 = b - Ax
	copy(b,r);
	A(x,Ap);
	x_plus_alpha_p(r,Ap,-1);

	// r0 = p0
	copy(r,p);

	for (size_t i = 0 ; i < 400 ; i++)
	{
		float alpha_c = alpha(A,r,p,Ap);

		x_plus_alpha_p(x,p,alpha_c);
		x_plus_alpha_p(rn,r,Ap,-alpha_c);

		float r_norm = norm(rn);

		std::cout << "It: " << i << "  " << r_norm << std::endl;

		if (r_norm < 0.1)
			return;

		float beta_c = beta(r,rn);
		x_plus_alpha_p(p,rn,p,beta_c);
		copy(rn,r);
	}
}



#endif /* EXAMPLE_NUMERICS_VORTEX_IN_CELL_CG_HPP_ */
