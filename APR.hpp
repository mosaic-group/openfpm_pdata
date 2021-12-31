//
#include <stdio.h>
#include <cstdlib>
#include "Operators/Vector/vector_dist_operators.hpp"
#include "OdeIntegrators/OdeIntegrators.hpp"
#include <iostream>
// #include "interpolation.hpp"
#include "Grid/grid_dist_id.hpp"
#include "Vector/vector_dist.hpp"
#include "Vector/vector_dist_subset.hpp"
#include "Decomposition/CartDecomposition.hpp"									   
#include "data_type/aggregate.hpp"
#include "DCPSE/DCPSE_op/DCPSE_op.hpp"
#include "DCPSE/SupportBuilder.hpp"
#include "APR.h"
#include "../alglib-cpp/src/interpolation.h"
#include "interpolation/mp4_kernel.hpp"
#include "interpolation/z_spline.hpp"
#include <boost/math/interpolators/vector_barycentric_rational.hpp>
#include <vector>
#include "APR.h"
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include "intr_2.h"
//
