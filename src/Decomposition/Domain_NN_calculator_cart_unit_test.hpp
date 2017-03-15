/*
 * Domain_NN_calculator_cart_unit_test.hpp
 *
 *  Created on: Mar 11, 2017
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_DOMAIN_NN_CALCULATOR_CART_UNIT_TEST_CPP_
#define SRC_DECOMPOSITION_DOMAIN_NN_CALCULATOR_CART_UNIT_TEST_CPP_

#include "VCluster/VCluster.hpp"
#include "Domain_NN_calculator_cart.hpp"
#include "NN/CellList/CellList.hpp"

BOOST_AUTO_TEST_SUITE (Domain_NN_cart_unit_test)

BOOST_AUTO_TEST_CASE( Domain_NN_cart_unit_test_usage )
{
	domain_nn_calculator_cart<3> dnn;

	Box<3,long int> pbox({20,30,40},{40,50,60});

	openfpm::vector<Box<3,size_t>> loc_box;
	loc_box.add(Box<3,size_t>({21,31,41},{30,40,50}));
	loc_box.add(Box<3,size_t>({31,41,51},{39,49,59}));

	size_t sz[3];
	sz[0] = pbox.getHigh(0) - pbox.getLow(0) + 1;
	sz[1] = pbox.getHigh(1) - pbox.getLow(1) + 1;
	sz[2] = pbox.getHigh(2) - pbox.getLow(2) + 1;

	grid_sm<3,void> gs2(sz);
	grid_key_dx<3> shift({5,5,5});

	dnn.setParameters(pbox);
	dnn.setNNParameters(loc_box,shift,gs2);

	// Linearized domain cell over pbox
	openfpm::vector_std<size_t> & cd = dnn.getDomainCells();
	openfpm::vector_std<size_t> & cdCSRdom = dnn.getCRSDomainCells();
	openfpm::vector_std<subsub_lin<3>> & cdCSRanom = dnn.getCRSAnomDomainCells();

	grid_sm<3,void> gs(sz);

	size_t tot_cl = 0;
	for (size_t i = 0 ; i < loc_box.size() ; i++)
		tot_cl += Box<3,size_t>(loc_box.get(i)).getVolumeKey();

	BOOST_REQUIRE_EQUAL(cd.size(),tot_cl);

	Box<3,size_t> box0 = loc_box.get(0);

	size_t dom_cell = (box0.getHigh(0) - box0.getLow(0)) *
					  (box0.getHigh(1) - box0.getLow(1)) *
					  (box0.getHigh(2) - box0.getLow(2) + 1);

	Box<3,size_t> box1 = loc_box.get(1);

	dom_cell += (box1.getHigh(0) - box1.getLow(0)) *
			  (box1.getHigh(1) - box1.getLow(1)) *
			  (box1.getHigh(2) - box1.getLow(2) + 1);

	BOOST_REQUIRE_EQUAL(cdCSRdom.size(),dom_cell);

	size_t anom_cell = 2*(box0.getHigh(0) - box0.getLow(0) + 1) * ((box0.getHigh(2) - box0.getLow(2)) + 1);
	anom_cell += 2*(box1.getHigh(0) - box1.getLow(0) + 1) * ((box1.getHigh(2) - box1.getLow(2)) + 1);

	anom_cell += 2*(box0.getHigh(1) - box0.getLow(1) + 1) * ((box0.getHigh(2) - box0.getLow(2)) + 1);
	anom_cell += 2*(box1.getHigh(1) - box1.getLow(1) + 1) * ((box1.getHigh(2) - box1.getLow(2)) + 1);

	BOOST_REQUIRE_EQUAL(cdCSRanom.size(),anom_cell);
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* SRC_DECOMPOSITION_DOMAIN_NN_CALCULATOR_CART_UNIT_TEST_CPP_ */
