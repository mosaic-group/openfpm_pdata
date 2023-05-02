//
// Created by jstark on 08.10.21
//

#ifndef DIFFUSION_HELPFUNCTIONSDIFFUSION_HPP
#define DIFFUSION_HELPFUNCTIONSDIFFUSION_HPP
#include "util/PathsAndFiles.hpp"

#include "level_set/redistancing_Sussman/HelpFunctionsForGrid.hpp"

/**@brief Get timestep that fulfills stability criterion nD-diffusion using a FTCS scheme
 *
 * @tparam grid_type Template type of input g_sparse.
 * @tparam T Template type of diffusion coefficient and timestep.
 * @param grid Input grid of type grid_type.
 * @param k_max Max. diffusion coefficient in the domain.
 * @return T type timestep.
 */
template <typename grid_type, typename T>
T diffusion_time_step(grid_type & grid, T k_max)
{
	T sum = 0;
	for(int d = 0; d < grid_type::dims; ++d)
	{
		sum += (T)1.0 / (T)(grid.spacing(d) * grid.spacing(d));
	}
	
	return 1.0 / ((T)4.0 * k_max * sum);
}

/**@brief Sigmoidal which can be shifted and stretched in order to get desired smooth inhomog. diffusion-coeff scalar
 * field
 *
 * @tparam T Template type of coordinates.
 * @param x T type abscissa axis of sigmoidal function.
 * @param x_shift T type shift of sigmoidal along abscissa.
 * @param x_stretch T type stretching of sigmoidal along acscissa - the smaller, the steeper!
 * @param y_min T type asymptotic lowest y-value.
 * @param y_max T type asymptotic maximal y-value.
 * @return Function value of sigmoidal for input value x.
 */
template <typename T>
T get_smooth_sigmoidal(T x, T x_shift, T x_stretch, T y_min, T y_max)
{
	return y_min + y_max / (1.0 + exp(- (x_shift + x_stretch * x)));
}

/**@brief Sums up property values over all particles in grid.
 *
 * @tparam Prop Index of property whose values will be summed up over all particles.
 * @tparam grid_type Template type of input particle vector.
 * @param grid Input particle vector. Can be a subset.
 * @return Sum of values stored in Prop.
 */
template <size_t Prop, typename grid_type>
auto sum_prop_over_grid(grid_type & grid)
{
	auto sum = 0.0;
	auto dom = grid.getDomainIterator();
	while(dom.isNext())
	{
		auto key = dom.get();
		sum += grid.template getProp<Prop>(key);
		++dom;
	}
	auto &v_cl = create_vcluster();
	v_cl.sum(sum);
	v_cl.execute();
	return sum;
}

/**@brief Writing out the total mass to a csv file.
 *
 * @tparam Conc Index of property containing the concentration.
 * @tparam T_mass Template type of initial mass.
 * @param grid Input particle vector_dist of type grid_type.
 * @param m_initial Initial total mass.
 * @param t Current time of the diffusion.
 * @param i Current iteration of the diffusion.
 * @param path_output std::string Path where csv file will be saved.
 * @param filename std::string Name of csv file. Default: "/total_mass.csv"
 * @Outfile 4 columns: Time, Total mass, Mass difference, Max mass
 */
template <size_t Conc, typename grid_type, typename T, typename T_mass>
void monitor_total_mass(grid_type & grid, const T_mass m_initial, const T p_volume, const T t, const size_t i,
                           const std::string & path_output, const std::string & filename="total_mass.csv")
{
	auto &v_cl = create_vcluster();
	
//	if (m_initial == 0)
//	{
//		if (v_cl.rank() == 0) std::cout << "m_initial is zero! Normalizing the total mass with the initial mass will "
//										   "not work. Aborting..." << std::endl;
//		abort();
//	}
	
	T m_total = sum_prop_over_grid<Conc>(grid) * p_volume;
	T m_diff = m_total - m_initial;
	T m_max = get_max_val<Conc>(grid) * p_volume;
	
	if (v_cl.rank() == 0)
	{
		std::string outpath = path_output + "/" + filename;
		create_file_if_not_exist(outpath);
		std::ofstream outfile;
		outfile.open(outpath, std::ios_base::app); // append instead of overwrite
		
		if (i == 0)
		{
			outfile   << "Time, Total mass, Mass difference, Max mass" << std::endl;
//			std::cout << "Time, Total mass, Total mass in-/decrease, Max mass" << std::endl;
		}
		
		outfile
				<< to_string_with_precision(t, 6) << ","
				<< to_string_with_precision(m_total, 6) << ","
				<< to_string_with_precision(m_diff, 6) << ","
				<< to_string_with_precision(m_max, 6)
				<< std::endl;
		outfile.close();
#if 0
		std::cout
				<< to_string_with_precision(t, 6) << ","
				<< to_string_with_precision(m_total, 6) << ","
				<< to_string_with_precision(m_diff, 6) << ","
				<< to_string_with_precision(m_max, 6)
				<< std::endl;
#endif
//		if (m_domain_normalized > 100 || m_domain_normalized < 0)
//		{
//			std::cout << "Mass increases or is negative, something went wrong. Aborting..." << std::endl;
//			abort();
//		}
	}
}

/**@brief Writing out the total concentration over time to a csv file.
 *
 * @tparam Conc Index of property containing the concentration.
 * @param grid Input particle vector_dist of type grid_type.
 * @param t Current time of the diffusion.
 * @param i Current iteration of the diffusion.
 * @param path_output std::string Path where csv file will be saved.
 * @param filename std::string Name of csv file. Default: "/normalized_total_mass.csv"
 * @Outfile 3 columns: timestep of writeout i, total concentration
 */
template <size_t Conc, typename grid_type, typename T>
void monitor_total_concentration(grid_type & grid, const T t, const size_t i,
                        const std::string & path_output, const std::string & filename="total_conc.csv")
{
	auto &v_cl = create_vcluster();
	T c_total = sum_prop_over_grid<Conc>(grid);

	if (v_cl.rank() == 0)
	{
		std::string outpath = path_output + "/" + filename;
		create_file_if_not_exist(outpath);
		std::ofstream outfile;
		outfile.open(outpath, std::ios_base::app); // append instead of overwrite

		if (i == 0)
		{
			outfile   << "Time, Total concentration" << std::endl;
		}

		outfile
				<< to_string_with_precision(t, 6) << ","
				<< to_string_with_precision(c_total, 6)
				<< std::endl;
		outfile.close();
	}
}

/**@brief Sets the emission term to 0 if concentration surpasses a threshold. This is to avoid accumulation of mass
 * at unconnected islands.
 *
 * @tparam Conc Index of property containing the concentration.
 * @tparam Emission Index of property containing the emission constant.
 * @tparam grid_type Template type of input particle vector.
 * @tparam key_vector_type Template type of particle keys.
 * @param grid Input particle vector. Can be a subset.
 * @param keys Keys for which necessity for adaptation is checked.
 * @param threshold Maximal concentration that is tolerated for keeping emission.
 */
template <size_t Conc, size_t Emission, typename grid_type, typename key_vector_type, typename T>
void adapt_emission(grid_type & grid, key_vector_type & keys, T threshold)
{
	for (int i = 0; i < keys.size(); ++i)
	{
		auto key = keys.get(i);
		if (grid.template getProp<Conc>(key) > threshold) grid.template getProp<Emission>(key) = 0.0;
	}
}

#endif //DIFFUSION_HELPFUNCTIONSDIFFUSION_HPP
