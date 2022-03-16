//
// Created by landfried on 14.03.22.
//

#ifndef OPENFPM_PDATA_INSTANCE_HPP
#define OPENFPM_PDATA_INSTANCE_HPP

#include "iostream"

template <typename ParticleMethodType, typename SimulationParametersType>
class Instance {

protected:

    void add() {
        std::cout << " add" << std::endl;
    }


public:

    virtual void freePlacement() {}


};


#endif //OPENFPM_PDATA_INSTANCE_HPP
