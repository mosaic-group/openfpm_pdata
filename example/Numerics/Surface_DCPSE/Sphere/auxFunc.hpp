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
// Created by foggia on 22-04-11
//


// 1) Function to compute analytically the SDF and the normal for the sphere in the whole domain
template <typename grid_type, int SDF, int normal>
void init_surfaceAndNormal(grid_type & grid,
			   const std::array<double,grid_type::dims> & center,
			   const double radius) {

  const int x{0};
  const int y{1};
  const int z{2};
  double dist, arg, theta, phi;
  std::array<double, grid_type::dims> coord;
  
  auto it = grid.getDomainIterator();
  while(it.isNext()) {
    auto key = it.get();

    dist = 0;
    for (int d = 0; d < grid_type::dims; ++d) {
      coord[d] = grid.getPos(key)[d];
      dist += (coord[d]-center[d])*(coord[d]-center[d]);
    }
    arg = (coord[x]-center[x]) * (coord[x]-center[x]) + (coord[y]-center[y]) * (coord[y]-center[y]);
    theta = std::atan2(std::sqrt(arg),(coord[z]-center[z]));
    phi = std::atan2((coord[y]-center[y]),(coord[x]-center[x]));
    
    grid.template get<SDF>(key) = std::sqrt(dist) - radius;
    grid.template get<normal>(key)[x] = std::sin(theta)*std::cos(phi);
    grid.template get<normal>(key)[y] = std::sin(theta)*std::sin(phi);
    grid.template get<normal>(key)[z] = std::cos(theta);
    
    ++it;
  }
}

// 2) Function to compute analytically the SDF and the normal for the sphere in the whole domain
template <typename grid_type, int SDF>
void init_surface(grid_type & grid,
		  const std::array<double,grid_type::dims> & center,
		  const double radius) {

  double dist;
  std::array<double, grid_type::dims> coord;
  
  auto it = grid.getDomainIterator();
  while(it.isNext()) {
    auto key = it.get();

    dist = 0;
    for (int d = 0; d < grid_type::dims; ++d) {
      coord[d] = grid.getPos(key)[d];
      dist += (coord[d]-center[d])*(coord[d]-center[d]);
    }
    
    grid.template get<SDF>(key) = std::sqrt(dist) - radius;    
    ++it;
  }
}

// 3a) Function to initialize the value of the quantity (spherical harmonic (4,0)) on the narrow band
// Y(4,0) = 3/16 * sqrt(1/pi) * (35*cos^4(theta) - 30*cos^2(theta) + 3)
template <typename grid_type, int qty, typename key_type>
void init_qty(grid_type & grid,
	      const std::array<double,grid_type::dims> & center,
	      std::vector<key_type> & nb_keys) {
  
  const int x{0};
  const int y{1};
  const int z{2};
  double arg, theta, cos2;
  const double pi{3.14159265358979323846};
  const double prefactor{3.0/16.0 * std::sqrt(1.0/pi)};
  std::array<double,grid.dims> coord;
  
  for (int i = 0; i < nb_keys.size(); ++i) {
    
    for (int d = 0; d < grid.dims; ++d)
      coord[d] = grid.getPos(nb_keys[i])[d];

    arg = (coord[x]-center[x]) * (coord[x]-center[x]) + (coord[y]-center[y]) * (coord[y]-center[y]);
    theta = std::atan2(std::sqrt(arg),(coord[z]-center[z]));
    cos2 = std::cos(theta)*std::cos(theta);
    
    grid.template get<qty>(nb_keys[i]) = prefactor *  (cos2*(35.0*cos2 - 30.0) + 3.0);
    
  }
}

// 3b) Function to initialize the value of the quantity (spherical harmonic (4,0)) on the whole grid
// Y(4,0) = 3/16 * sqrt(1/pi) * (35*cos^4(theta) - 30*cos^2(theta) + 3)
template <typename grid_type, int qty>
void init_qty(grid_type & grid,
	      const std::array<double,grid_type::dims> & center) {
  
  const int x{0};
  const int y{1};
  const int z{2};
  double arg, theta, cos2;
  const double pi{3.14159265358979323846};
  const double prefactor{3.0/16.0 * std::sqrt(1.0/pi)};
  std::array<double,grid.dims> coord;

  auto it = grid.getDomainIterator();
  while(it.isNext()) {
    auto key = it.get();

    for (int d = 0; d < grid.dims; ++d)
      coord[d] = grid.getPos(key)[d];
    
    arg = (coord[x]-center[x]) * (coord[x]-center[x]) + (coord[y]-center[y]) * (coord[y]-center[y]);
    theta = std::atan2(std::sqrt(arg),(coord[z]-center[z]));
    cos2 = std::cos(theta)*std::cos(theta);
    
    grid.template get<qty>(key) = prefactor *  (cos2*(35.0*cos2 - 30.0) + 3.0);
    ++it;
  }
}
// 4) Function to compute the analytical solution of the LB operator in the whole domain
template <typename grid_type, int qty, int sol>
void init_analytSol(grid_type & grid,
		    const std::array<double,grid_type::dims> & center,
		    const double radius) {

  const int x{0};
  const int y{1};
  const int z{2};
  double arg, theta, cos2;
  const double pi{3.14159265358979323846};
  const double prefactor{3.0/16.0 * std::sqrt(1.0/pi)};
  std::array<double,grid.dims> coord;

  auto it = grid.getDomainIterator();
  while(it.isNext()) {
    auto key = it.get();
    
    for (int d = 0; d < grid.dims; ++d)
      coord[d] = grid.getPos(key)[d];
    
    arg = (coord[x]-center[x]) * (coord[x]-center[x]) + (coord[y]-center[y]) * (coord[y]-center[y]);
    theta = std::atan2(std::sqrt(arg),(coord[z]-center[z]));
    cos2 = std::cos(theta)*std::cos(theta);
    
    grid.template get<sol>(key) = -20.0 * prefactor * (cos2*(35.0*cos2 - 30.0) + 3.0);   
    ++it;
  }
}

// 5) Functions to get the points on a NB around the surface (from jstark - Sussman code)
bool within_narrow_band(double phi,
			double b_low,
			double b_up) {
  return (phi >= b_low && phi <= b_up);
}

template<typename grid_type, size_t prop, typename key_type>
void get_NB_indices(grid_type & grid,
                    double thickness, // NOT grid points
                    std::vector<key_type> & nb_keys) {

  double b_low = -thickness/2.0;
  double b_up  =  thickness/2.0;

  auto it{grid.getDomainIterator()};
  while (it.isNext()) {

    auto key{it.get()};
    if (within_narrow_band(grid.template get<prop>(key),b_low,b_up))
      nb_keys.push_back(key);
    ++it;
  }
}

// 6) Functions to compute the L2 and Linf norms (from jstark)
struct L_norms {
  double l2;
  double linf;
};

template<typename grid_type, int PropNumeric, int PropAnalytic, int Error, typename key_type>
void get_absolute_error(grid_type & grid,
			std::vector<key_type> & nb_keys) {
  for (int i = 0; i < nb_keys.size(); ++i)
    grid.template getProp<Error>(nb_keys[i]) = std::fabs(grid.template getProp<PropAnalytic>(nb_keys[i]) - (grid.template getProp<PropNumeric>(nb_keys[i])));
}

template<typename grid_type, int error, typename key_type>
L_norms get_l_norms_NB(grid_type & grid,
		       std::vector<key_type> & nb_keys) {

  double maxError{0};
  double sumErrorSq{0};
  size_t num_nb{nb_keys.size()};

  for (int i = 0; i < nb_keys.size(); ++i) {
    
    sumErrorSq += grid.template getProp<error>(nb_keys[i]) * grid.template getProp<error>(nb_keys[i]);
    if (grid.template getProp<error>(nb_keys[i]) > maxError)
      maxError = grid.template getProp<error>(nb_keys[i]); // update maxError
  }

  // Communicate between processes
  auto &v_cl = create_vcluster();
  v_cl.sum(num_nb);
  v_cl.sum(sumErrorSq);
  v_cl.max(maxError);
  v_cl.execute();
  
  double l2{std::sqrt(sumErrorSq / (double)num_nb)};
  // std::cout << "sum of num_nb: " << (double)num_nb << "\n";
  double linf{maxError};
  return {l2, linf};
}

template<typename T>
void write_lnorms_to_file(T N, L_norms l_norms, std::string filename, std::string path_output) {
  auto &v_cl = create_vcluster();
  if (v_cl.rank() == 0) {
    std::string path_output_lnorm{path_output + "/" + filename + ".csv"};
    create_file_if_not_exist(path_output_lnorm);
    
    std::ofstream l_out;
    l_out.open(path_output_lnorm, std::ios::app); // append instead of overwrite
    l_out << N << ',' << std::setprecision(6) << std::scientific << l_norms.l2 << ',' << l_norms.linf
	  << std::endl;
    l_out.close();
  }
}

// Function to set particle property to zero
template<typename vector_type, int prop>
void set_prop2zero(vector_type & vec) {
  auto it = vec.getDomainIterator();
  
  while (it.isNext()) {
    auto key = it.get();
    vec.template getProp<prop>(key) = 0.0;    
    ++it;
  }
}


// 7) Function to initialize the projection matrix of the surface
// template<typename grid_type, int normal, int projMat_prop>
// void init_projMat(grid_type & grid) {

//   double ni, nj;
//   auto it{grid.getDomainIterator()};
//   while (it.isNext()) {
//     auto key{it.get()};
//     for (int i = 0; i < DIM; ++i) {
//       ni = grid.template get<normal>(key)[i];
//       for (int j = 0; j < DIM; ++j) {
// 	nj = grid.template get<normal>(key)[j];
// 	grid.template get<projMat_prop>(key)[i][j] = (i==j)*(1-ni*nj) - !(i==j)*(ni*nj);
//       }
//     }
//     ++it;
//   }
// }


// 8) Function to compute the surface gradient of a quantity
// template<typename grid_type, int QTY, int GRAD, int SURFGRAD, int PROJMAT, typename DX, typename DY, typename DZ, typename key_type>
// void SGrad(grid_type & grid,
// 	   DX & Dx,
// 	   DY & Dy,
// 	   DZ & Dz) {

//   auto qty{FD::getV<QTY>(grid)};
//   auto grad_qty{FD::getV<GRAD>(grid)};
  
//   grid.template ghost_get<QTY,GRAD,PROJMAT>();

//   grad_qty[0] = Dx(qty);
//   grad_qty[1] = Dy(qty);
//   grad_qty[2] = Dz(qty);

//   // Projection
//   auto it{grid.getDomainIterator()};
//   while (it.isNext()) {
    
//     auto key{it.get()};
    
//     for (int l = 0; l < DIM; ++l)
//       for (int k = 0; k < DIM; ++k)
// 	for (int t = 0; t < DIM; ++t)
// 	  for (int h = 0; h < DIM; ++h)
// 	    grid.template get<SURFGRAD>(key)[l][k] += grid.template get<PROJMAT>(key)[l][t] * grid.template get<PROJMAT>(key)[k][h] * grid.template get<EUCGRAD>(key)[t][h]; 
//     ++it;
//   }
  
// }
