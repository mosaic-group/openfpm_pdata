//
// Created by landfried on 03.02.22.
//

#include <functional>

#ifndef OPENFPM_PDATA_GHOST_HPP
#define OPENFPM_PDATA_GHOST_HPP


template <int ... prp, typename ParticleDataType>
void ghost_get_N (ParticleDataType& particleData, std::integer_sequence<int, prp...>)
{
    particleData.getOpenFPMContainer().template ghost_get<prp...>();
}




#endif //OPENFPM_PDATA_GHOST_HPP
