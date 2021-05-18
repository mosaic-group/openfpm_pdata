
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Space/Shape/Box.hpp"
#include "Decomposition/Domain_icells_cart.hpp"
#include "VCluster/VCluster.hpp"
#include "VTKWriter/VTKWriter.hpp"
#include "Grid/iterators/grid_skin_iterator.hpp"

BOOST_AUTO_TEST_SUITE( domain_icells_cart )

BOOST_AUTO_TEST_CASE( domain_icells_use )
{
#if 0

	domain_icell_calculator<3,float,memory_traits_inte,CudaMemory> dcc;

	openfpm::vector_gpu<SpaceBox<3,float>> domain_proc;

	Box<3,float> box1({0.0,0.0,0.0},{0.5,0.5,0.5});
	Box<3,float> box2({0.5,0.25,0.25},{0.75,0.5,0.5});

	domain_proc.add(box1);
	domain_proc.add(box2);

	openfpm::vector_gpu<Box<3,float>> ie_ghost;

	// first sub ghost area

	Box<3,float> g1({0.0,0.0,0.0},{0.1,0.5,0.5});
	Box<3,float> g2({0.0,0.0,0.4},{0.5,0.5,0.5});
	Box<3,float> g3({0.0,0.0,0.0},{0.5,0.1,0.5});
	Box<3,float> g4({0.0,0.4,0.0},{0.5,0.5,0.5});
	Box<3,float> g5({0.0,0.0,0.0},{0.5,0.5,0.1});
	Box<3,float> g6_1({0.4,0.0,0.0},{0.5,0.25,0.5});
	Box<3,float> g6_2({0.4,0.0,0.0},{0.5,0.5,0.25});

	// Second sub ghost area

	Box<3,float> g7({0.5,0.25,0.25},{0.6,0.5,0.5});
	Box<3,float> g8({0.65,0.25,0.25},{0.75,0.5,0.5});
	Box<3,float> g9({0.5,0.25,0.25},{0.75,0.35,0.25});
	Box<3,float> g10({0.5,0.4,0.25},{0.75,0.5,0.5});
	Box<3,float> g11({0.5,0.25,0.25},{0.75,0.5,0.35});
	Box<3,float> g12({0.5,0.25,0.4},{0.75,0.5,0.5});

	ie_ghost.add(g1);
	ie_ghost.add(g2);
	ie_ghost.add(g3);
	ie_ghost.add(g4);
	ie_ghost.add(g5);
	ie_ghost.add(g6_1);
	ie_ghost.add(g6_2);

	ie_ghost.add(g7);
	ie_ghost.add(g8);
	ie_ghost.add(g9);
	ie_ghost.add(g10);
	ie_ghost.add(g11);
	ie_ghost.add(g12);

	Box<3,float> pbox({0.0,0.0,0.0},{0.75,0.5,0.5});

	Ghost<3,float> g(0.1);

	auto & v_cl = create_vcluster<CudaMemory>();

	dcc.CalculateInternalCells(v_cl,ie_ghost,domain_proc,pbox,0.1,g);

	// Check

	auto & ic = dcc.getIcells();
	auto & dc = dcc.getDcells();

	ic.template deviceToHost<0>();
	dc.template deviceToHost<0>();

	openfpm::vector<unsigned int> icheck;
	openfpm::vector<unsigned int> dcheck;

	grid_sm<3,void> gr = dcc.getGrid();

	Box<3,size_t> b1({2,2,2},{6,6,6});

	size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	grid_skin_iterator_bc<3> it_skin(gr,b1,b1,bc);

	while (it_skin.isNext())
	{
		auto p = it_skin.get();

		if (p.get(0) == 6 && p.get(1) == 5 && p.get(2) == 5)
		{}
		else
		{icheck.add(gr.LinId(p));}

		++it_skin;
	}

	Box<3,size_t> b3({7,4,4},{9,6,6});
	grid_key_dx_iterator_sub<3> it_sub(gr,b3.getKP1(),b3.getKP2());

	while (it_sub.isNext())
	{
		auto p = it_sub.get();

		icheck.add(gr.LinId(p));

		++it_sub;
	}

	// dcheck has only the internal part of the first domain + 6,5,5

	Box<3,size_t> b1i({3,3,3},{5,5,5});
	grid_key_dx_iterator_sub<3> it_sub2(gr,b1i.getKP1(),b1i.getKP2());

	while (it_sub2.isNext())
	{
		auto p = it_sub2.get();

		dcheck.add(gr.LinId(p));

		++it_sub2;
	}

	grid_key_dx<3> kk({6,5,5});
	dcheck.add(gr.LinId(kk));

	icheck.sort();
	dcheck.sort();

	for (size_t i = 0 ; i < icheck.size() ; i++)
	{BOOST_REQUIRE_EQUAL(icheck.template get<0>(i),ic.template get<0>(i));}

	for (size_t i = 0 ; i < dcheck.size() ; i++)
	{BOOST_REQUIRE_EQUAL(dcheck.template get<0>(i),dc.template get<0>(i));}

	#endif
}


BOOST_AUTO_TEST_SUITE_END()
