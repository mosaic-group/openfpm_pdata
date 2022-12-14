//
// Created by landfried on 11.11.22.
//

#ifndef OPENFPM_PDATA_ALIAS_HPP
#define OPENFPM_PDATA_ALIAS_HPP



#define PARTICLE(property_arg) particle.template property<property_arg>()
#define NEIGHBOR(property_arg) neighbor.template property<property_arg>()
#define NEW_PARTICLE(property_arg) particle.template new_property<property_arg>()



#endif //OPENFPM_PDATA_ALIAS_HPP
