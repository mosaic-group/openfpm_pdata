/*
 * DLB.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Antonio Leo
 */

#ifndef SRC_DECOMPOSITION_DLB_HPP_
#define SRC_DECOMPOSITION_DLB_HPP_

//! Time structure for statistical purposes
typedef struct
{
	//! starting time of the simulation (0)
	size_t simulationStartTime = 0;

	//! End iteration of the simulation
	size_t simulationEndTime;

	//! integration time
	double timeStep = 0.1;

	//! Interval between teo rebalance

	//! Start time
	size_t iterationStartTime;

	//! End time
	size_t iterationEndTime;
} Times;

/*! Class that implements the two heuristics to determine when a re-balance of the distribution is needed.
 *
 *  Used heuristics are: SAR and Un-balance Threshold (Default)\n
 *
 *  To chose the heuristic use the method setHeuristic(Heuristic)
 *
 *  In the SAR heuristic the following formula is applied:\n
 *  \f$W_{n} = \frac{\sum_{j=1}^{n} (T_{max}(j) - T_{avg}(j)) + C} {n}\f$
 *
 *  \f$T_{max}(j)\f$ – wall-clock time of bottleneck process in time step j\n
 *  \f$T_{avg}(j)\f$ – average wall-clock time for time step j over all processes\n
 *  \f$C\f$ – cost of re-decomposing the problem\n
 *  \f$n\f$ – number of time steps since last re-decomposition\n
 *  \n
 *  For small n, load balance is good and W decreases since C is amortized over an increasing number of time steps.
 *  As the accumulated idle time starts to dominate, W starts to rise. At this point, C has been fully amortized.
 *  Re-decompose when \f$W_{n} > W_{n-1}\f$\n
 *
 *  In the Un-balance Threshold heuristic the re-balance is triggered when the un-balance level exceeds a certain level.
 *  Levels can be chosen in the ThresholdLevel type.
 */
class DLB
{
public:

	//! Type of DLB heuristics
	enum Heuristic
	{
		SAR_HEURISTIC, UNBALANCE_THRLD
	};

	//! Level of un-balance needed to trigger the re-balance
	enum ThresholdLevel
	{
		THRLD_LOW = 5, THRLD_MEDIUM = 7, THRLD_HIGH = 10
	};

private:

	//! Runtime virtual cluster machine
	Vcluster & v_cl;

	//! Structure that will contain all the timings
	Times timeInfo;

	//! Wn for SAR heuristic
	float w_n = -1;

	//! Computation cost for SAR heuristic
	float c_c = 5;

	//! Number of time-steps since the previous DLB
	size_t n_ts = 1;

	//! Idle time accumulated so far, needed for SAR heuristic
	float i_time = 0;

	//! Vector to collect all timings
	openfpm::vector<long> times;

	//! Type of the heuristic to use
	Heuristic heuristic = UNBALANCE_THRLD;

	//! Un-balance value
	float unbalance = -1;

	//! Threshold value
	ThresholdLevel thl = THRLD_MEDIUM;

	/*! \brief Function that gather times informations and decides if a rebalance is needed it uses the SAR heuristic
	 *
	 * \return true if re-balance is needed
	 *
	 */
	inline bool SAR()
	{
		long t = timeInfo.iterationEndTime - timeInfo.iterationStartTime;
		float t_max = t, t_avg = t;

		// Exchange time informations through processors
		v_cl.max(t_max);
		v_cl.sum(t_avg);
		v_cl.execute();

		t_avg /= v_cl.getProcessingUnits();

		// add idle time to vector
		i_time += t_max - t_avg;

		// Compute Wn
		float nw_n = (i_time + c_c) / n_ts;

		if (w_n == -1)
			w_n = nw_n;

		if (nw_n > w_n)
		{
			i_time = 0;
			n_ts = 1;
			w_n = nw_n;
			return true;
		}
		else
		{
			++n_ts;
			w_n = nw_n;
			return false;
		}
	}

	/*! \brief Check if the un-balance has exceeded the threshold
	 *
	 * \return true if re-balance is needed, false otherwise
	 */
	bool unbalanceThreshold()
	{
		if (unbalance == -1)
		{
			std::cerr << "Error: Un-balance value must be set before checking DLB.";
			return false;
		}

		if (unbalance > thl)
		{
			return true;
		}

		return false;
	}

public:

	/*! \brief Constructor for DLB class
	 *
	 * \param v_cl virtual cluster object
	 */
	DLB(Vcluster & v_cl) :
			v_cl(v_cl)
	{
	}

	/*! \brief Set the heuristic to use (default: un-balance threshold)
	 *
	 * \param h
	 */
	void setHeurisitc(Heuristic h)
	{
		heuristic = h;
	}

	/*! \brief Get the heuristic
	 *
	 * Indicate which heuristic model is used to calculate when a rebalance
	 * is needed
	 *
	 * \return the Heuristic used by DLB
	 *
	 */
	Heuristic getHeurisitc()
	{
		return heuristic;
	}

	/*! \brief check if a re-balance is needed using the selected heuristic
	 *
	 * \return true if the rebalance is needed
	 *
	 */
	bool rebalanceNeeded()
	{
		if (heuristic == SAR_HEURISTIC)
		{
			return SAR();
		}
		else
		{
			return unbalanceThreshold();
		}
	}

	/*! \brief Set start time for the simulation
	 *
	 * \param t time when the whole simulation starts
	 */
	void setSimulationStartTime(size_t t)
	{
		timeInfo.simulationStartTime = t;
	}

	/*! \brief Get start time for the simulation
	 *
	 * \return the start point of the simulation
	 *
	 */
	size_t getSimulationStartTime()
	{
		return timeInfo.simulationStartTime;
	}

	/*! \brief Set end time for the simulation
	 *
	 * \param t time when the whole simulation ends
	 */
	void setSimulationEndTime(size_t t)
	{
		timeInfo.simulationEndTime = t;
	}

	/*! \brief Get end time for the simulation
	 *
	 * \return the end time of the simulation
	 *
	 */
	size_t getSimulationEndTime()
	{
		return timeInfo.simulationEndTime;
	}

	/*! \brief Set start time for the single iteration
	 *
	 */
	void startIteration()
	{
		timeInfo.iterationStartTime = clock();
	}

	/*! \brief Set start time for the single iteration
	 *
	 * \param t time when the one iteration starts
	 */
	void startIteration(size_t t)
	{
		timeInfo.iterationStartTime = t;
	}

	/*! \brief Set end time for the single iteration
	 *
	 * \param time when one iteration is completed
	 *
	 */
	void endIteration()
	{
		timeInfo.iterationEndTime = clock();
	}

	/*! \brief Set the end time when the previous rebalance has been performed
	 *
	 * \param t time when one iteration ends
	 */
	void endIteration(size_t t)
	{
		timeInfo.iterationEndTime = t;
	}

	/*! \brief Set delta time step for one iteration (Computation time)
	 *
	 * \param t timestep
	 */
	void setTimeStep(double t)
	{
		timeInfo.timeStep = t;
	}

	/*! \brief Set time step for the single iteration
	 *
	 * \param computation value of the computation cost (default: 5)
	 */
	void setComputationCost(size_t computation)
	{
		c_c = computation;
	}

	/*! \brief Get how many time-steps have passed since the last re-balancing
	 *
	 * \return number of timesteos
	 *
	 */
	size_t getNTimeStepSinceDLB()
	{
		return n_ts;
	}

	/*! \brief Set un-balance value
	 *
	 * \param u unbalance
	 */
	void setUnbalance(float u)
	{
		unbalance = u;
	}

	/*! \brief threshold of umbalance to start a rebalance
	 *
	 * \param t threshold level
	 */
	void setThresholdLevel(ThresholdLevel t)
	{
		thl = t;
	}

};

#endif /* SRC_DECOMPOSITION_DLB_HPP_ */
