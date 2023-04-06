// ------------------------------------------------------------------------------------------
//    ActiveGelsOnDeformableSurfaces
//    Copyright (C) 2020 ---- foggia
//
//    This file is part of the ActiveGelsOnDeformableSurfaces software.
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ------------------------------------------------------------------------------------------
//
// Created by foggia on 22-04-12
//

// Include standard library files
#include <iostream>

// Include OpenFPM files
#include "Grid/grid_dist_id.hpp"
#include "Vector/vector_dist.hpp"
#include "FiniteDifference/FD_op.hpp"
#include "level_set/closest_point/closest_point.hpp"
#include <util/PathsAndFiles.hpp>
#include "interpolation/interpolation.hpp"
#include "interpolation/mp4_kernel.hpp"

// Include own files
#include "auxFunc.hpp"
#include "spatial_discretization/finite_difference/derivative.hpp"
#include "interp_kernel/m4_kernel.hpp"

// ------------------------------------------------------------------------------------------
// Main function
// ------------------------------------------------------------------------------------------

int main(int argc, char * argv[]) {

  openfpm_init(&argc,&argv);
  Vcluster<> & v_cl = create_vcluster();

  timer ttotal, tcp, tExtQ, tExtRHS, tLap, tAll;
  const size_t mpiSize{v_cl.size()};
  
  const int DIM{3};
  const int SDF{0};
  const int CP{1};
  const int QTY{2};
  const int ANALYTLB{3};
  const int NUMLB{4};
  const int ABSERR{5};
  const int TMP_EXT_SCAL{6};

  const int NUMLB_PART{0};
  const int ANALYTLB_PART{1};
  const int ABSERR_PART{2};
  
  // Writing directory
  std::string cwd{boost::filesystem::current_path().string()};
  std::string output_dir{cwd + "/output_sphereR1_CartLap"};
  create_directory_if_not_exist(output_dir);
  std::string g_output;
  std::ofstream norms_file, normsParticles_file, timesCPU_file, timesWCT_file;
  if (v_cl.rank() == 0) {
    norms_file.open(output_dir + "/norms_EucLap_mpiSize_" + std::to_string(mpiSize) + ".csv", std::ios_base::app);
    normsParticles_file.open(output_dir + "/normsParticles_EucLap_mpiSize_" + std::to_string(mpiSize) + ".csv", std::ios_base::app);
    timesCPU_file.open(output_dir + "/timesCPU_EucLap_mpiSize_" + std::to_string(mpiSize) + ".csv", std::ios_base::app);
    timesWCT_file.open(output_dir + "/timesWCT_EucLap_mpiSize_" + std::to_string(mpiSize) + ".csv", std::ios_base::app);
    if (norms_file.is_open() == false or normsParticles_file.is_open() == false or timesCPU_file.is_open() == false or timesWCT_file.is_open() == false) {
      std::cout << "ERROR: File(s) is(are) not open.\n";
      return 1;
    }
    norms_file << "#N   L2   Linf" << std::endl;
    normsParticles_file << "#N   L2   Linf" << std::endl;
    timesCPU_file << "#N   tcp   tExtQ   tAll   tExtRHS   tLap   ttotal" << std::endl;
    timesWCT_file << "#N   tcp   tExtQ   tAll   tExtRHS   tLap   ttotal" << std::endl;
  }
  
  // Surface
  const double radius{1.0};
  const std::array<double,DIM> center{0.0,0.0,0.0};
  const int nb_width{10};
  const int nb_half_width{nb_width/2};
  const int POLYORD{3};
  const int POLYORD_EXT{3};
  
  // Domain
  // double boxP1{-1.3}, boxP2{1.3};
  double boxP1{-2.0}, boxP2{2.0};
  double boxSize{boxP2-boxP1};
  Box<DIM,double> domain{{boxP1,boxP1,boxP1},{boxP2,boxP2,boxP2}};
  Ghost<DIM,long int> ghost{20};
  periodicity<DIM> bc{{NON_PERIODIC,NON_PERIODIC,NON_PERIODIC}};
  typedef aggregate<double,double[DIM],double,double,double,double,double> prop;
  typedef grid_dist_id<DIM,double,prop> grid_type;
  openfpm::vector<std::string> GpropNames{"phiSDF","cp","qty","analytLB","numLB","absErr","tmp_ext_scal"};

  // Particles on surface
  typedef aggregate<double,double,double> prop_part;
  openfpm::vector<std::string> PpropNames{"numLB","analytLB","absErr"};
  typedef vector_dist<DIM,double,prop_part> vector_type;

  std::vector<grid_dist_key_dx<grid_type::dims,grid_key_dx<grid_type::dims>>> nb_keys, nb_keys_error;

  int n[7] = {73,87,110,140,178,211,250};
  int n_part[7] = {3600,6000,12000,24000,48000,80000,128000};
  // int n[7] = {23,31,41,49,59,67,77};
  // int n[2] = {41,53};
  // int n[6] = {21,41,81,161,321,641};
  
  // Loop over different grid sizes
  for (int i = 6; i < 7; ++i) {

    ttotal.reset();
    tcp.reset();
    tExtQ.reset();
    tExtRHS.reset();
    tAll.reset();
    tLap.reset();

    g_output = {output_dir + "/sphereR1_SurfLap_" + std::to_string(n[i]) + "_mpiSize_" + std::to_string(mpiSize)};

    // Grid creation
    size_t sz[DIM] = {n[i],n[i],n[i]};
    grid_type grid{sz,domain,ghost,bc};
    grid.setPropNames(GpropNames);

    // Particle creation
    vector_type particlesSurf{grid.getDecomposition(),0};
    particlesSurf.setPropNames(PpropNames);

    std::array<double,DIM> coord;
    // const int n_part{4*n[i]};
    if (v_cl.rank() == 0) {

      // Created using the Fibonacci sphere algorithm
      const double pi{3.14159265358979323846};
      const double golden_ang = pi*(3.0 - std::sqrt(5.0));
      const double prefactor{3.0/16.0 * std::sqrt(1.0/pi)};
      double rad, theta, arg, thetaB, cos2;
    
      for (int j = 0; j < n_part[i]; ++j) {
      
	coord[1] = 1.0 - 2.0*(j/double(n_part[i]-1));
	rad = std::sqrt(1.0 - (coord[1]-center[1])*(coord[1]-center[1]));
	theta = golden_ang * j;
	coord[0] = std::cos(theta) * rad;
	coord[2] = std::sin(theta) * rad;

	arg = (coord[0]-center[0]) * (coord[0]-center[0]) + (coord[1]-center[1]) * (coord[1]-center[1]);
	thetaB = std::atan2(std::sqrt(arg),(coord[2]-center[2]));
	cos2 = std::cos(thetaB)*std::cos(thetaB);
      
	particlesSurf.add();
	particlesSurf.getLastPos()[0] = coord[0];
	particlesSurf.getLastPos()[1] = coord[1];
	particlesSurf.getLastPos()[2] = coord[2];
	particlesSurf.getLastProp<ANALYTLB_PART>() = -20.0 * prefactor * (cos2*(35.0*cos2 - 30.0) + 3.0); ;
      }
    }
    particlesSurf.map();

    if (v_cl.rank() == 0)
      std::cout << "n: " << n[i]
		<< " dx: " << grid.spacing(0)
		<< std::endl;

    // SDF and normal
    init_surface<grid_type,SDF>(grid,center,radius);
    grid.template ghost_get<SDF>();

    // NB indices
    nb_keys.clear();
    nb_keys_error.clear();
    get_NB_indices<grid_type,SDF>(grid,7*grid.spacing(0),nb_keys);
    get_NB_indices<grid_type,SDF>(grid,2*grid.spacing(0),nb_keys_error);
    
    // Analytical solution
    init_analytSol<grid_type,QTY,ANALYTLB>(grid,center,radius);

    // grid.write_frame(g_output,0,VTK_WRITER);
    // openfpm_finalize();
    // return 0;
    // -------------------------------------------------------------
    ttotal.start();
    
    // 1) CP representation
    tcp.start();
    estimateClosestPoint<SDF,CP,POLYORD>(grid,4*grid.spacing(0));
    grid.template ghost_get<CP>();
    tcp.stop();

    tAll.start();
    // 2) Quantity + extension
    init_qty<grid_type,QTY>(grid,center,nb_keys);
    grid.template ghost_get<QTY>();
    tExtQ.start();
    extendLSField<SDF,CP,QTY,TMP_EXT_SCAL,3>(grid,4*grid.spacing(0));
    grid.template ghost_get<QTY>();
    tExtQ.stop();
    // init_qty<grid_type,QTY>(grid,center);
    // grid.template ghost_get<QTY>();

    // 3) Euclidean Laplacian
    tLap.start();
    auto it = grid.getDomainIterator();
    while(it.isNext()) {
      auto key = it.get();
      grid.template get<NUMLB>(key) = second_order_deriv2<grid_type,QTY>(grid,key,0) + second_order_deriv2<grid_type,QTY>(grid,key,1) + second_order_deriv2<grid_type,QTY>(grid,key,2);
      ++it;
    }
    tLap.stop();
    // FD::Derivative<0,2,2,FD::CENTRAL> Dxx;
    // FD::Derivative<1,2,2,FD::CENTRAL> Dyy;
    // FD::Derivative<2,2,2,FD::CENTRAL> Dzz;
    // auto qty{FD::getV<QTY>(grid)};
    // auto numLB{FD::getV<NUMLB>(grid)};
    // tLap.start();
    // numLB = Dxx(qty) + Dyy(qty) + Dzz(qty);
    // tLap.stop();

    // 4) Extend Laplacian
    tExtRHS.start();
    grid.template ghost_get<NUMLB>();
    extendLSField<SDF,CP,NUMLB,TMP_EXT_SCAL,3>(grid,4*grid.spacing(0));
    grid.template ghost_get<NUMLB>();
    tExtRHS.stop();

    tAll.stop();
    
    ttotal.stop();
    // -------------------------------------------------------------
    
    std::cout << "n: " << n[i] << " h: " << grid.spacing(0) << std::endl;

    // Error on grid
    L_norms abs_norm_nb;
    get_absolute_error<grid_type,NUMLB,ANALYTLB,ABSERR>(grid,nb_keys_error);
    abs_norm_nb = get_l_norms_NB<grid_type,ABSERR>(grid,nb_keys_error);
    write_lnorms_to_file(n[i],abs_norm_nb,"norms_EucLap_mpiSize_" + std::to_string(mpiSize),output_dir);

    // Interpolate to particles
    typedef interpolate<vector_type,grid_type,mp4_kernel<double>> interp_type;
    interp_type interpSurf(particlesSurf,grid);
    set_prop2zero<vector_type,0>(particlesSurf);
    interpSurf.template m2p<NUMLB,NUMLB_PART>(grid,particlesSurf);
    
    // Error on particles
    double maxError{0}, sumErrorSq{0};
    auto pit{particlesSurf.getDomainIterator()};
    while(pit.isNext()) {
      auto key{pit.get()};

      particlesSurf.template getProp<ABSERR_PART>(key) = std::fabs(particlesSurf.template getProp<ANALYTLB_PART>(key) - particlesSurf.template getProp<NUMLB_PART>(key));
      sumErrorSq += particlesSurf.template getProp<ABSERR_PART>(key)*particlesSurf.template getProp<ABSERR_PART>(key);
      maxError = std::max(maxError,particlesSurf.template getProp<ABSERR_PART>(key));
      
      ++pit;
    }
    v_cl.sum(sumErrorSq);
    v_cl.max(maxError);
    v_cl.execute();
    abs_norm_nb.l2 = {std::sqrt(sumErrorSq / (double)n_part[i])};
    abs_norm_nb.linf = maxError;
    write_lnorms_to_file(n_part[i],abs_norm_nb,"normsParticles_EucLap_mpiSize_" + std::to_string(mpiSize),output_dir);
    
    // Time computation
    double tcp_cpu{tcp.getcputime()}, tExtQ_cpu{tExtQ.getcputime()}, tExtRHS_cpu{tExtRHS.getcputime()}, tLap_cpu{tLap.getcputime()}, tAll_cpu{tAll.getcputime()}, ttotal_cpu{ttotal.getcputime()};
    double tcp_wct{tcp.getwct()}, tExtQ_wct{tExtQ.getwct()}, tExtRHS_wct{tExtRHS.getwct()}, tLap_wct{tLap.getwct()}, tAll_wct{tAll.getwct()}, ttotal_wct{ttotal.getwct()};
    v_cl.sum(tcp_cpu);
    v_cl.sum(tExtQ_cpu);
    v_cl.sum(tExtRHS_cpu);
    v_cl.sum(tLap_cpu);
    v_cl.sum(tAll_cpu);
    v_cl.sum(tcp_wct);
    v_cl.sum(tExtQ_wct);
    v_cl.sum(tExtRHS_wct);
    v_cl.sum(tLap_wct);
    v_cl.sum(tAll_wct);
    v_cl.sum(ttotal_cpu);
    v_cl.sum(ttotal_wct);
    
    v_cl.execute();
    // grid.write_frame(g_output,0,VTK_WRITER);

    if (v_cl.rank() == 0) {
      timesCPU_file << n[i] << "," << std::setprecision(6) << std::scientific << tcp_cpu/v_cl.size() << "," <<  tExtQ_cpu/v_cl.size() << "," <<  tAll_cpu/v_cl.size() << "," <<  tExtRHS_cpu/v_cl.size() << "," <<  tLap_cpu/v_cl.size() << "," <<  ttotal_cpu/v_cl.size() << std::endl;
      timesWCT_file << n[i] << "," << std::setprecision(6) << std::scientific << tcp_wct/v_cl.size() << "," <<  tExtQ_wct/v_cl.size() << "," <<  tAll_wct/v_cl.size() << "," <<  tExtRHS_wct/v_cl.size() << "," <<  tLap_wct/v_cl.size() << "," <<  ttotal_wct/v_cl.size() << std::endl;
    }
  }
  norms_file.close();

  openfpm_finalize();
  return 0;
}
