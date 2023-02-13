//
// Created by landfried on 11.11.22.
//

#ifndef OPENFPM_PDATA_ALIAS_HPP
#define OPENFPM_PDATA_ALIAS_HPP



#define PARTICLE(property_arg) particle.template property_point<property_arg>()
#define NEIGHBOR(property_arg) neighbor.template property_point<property_arg>()
#define PARTICLE_WRITE(property_arg) particle.template property<property_arg>()
#define NEIGHBOR_WRITE(property_arg) neighbor.template property<property_arg>()


#define PARTICLE_POSITION particle.position_point()
#define NEIGHBOR_POSITION neighbor.position_point()
#define PARTICLE_POSITION_WRITE particle.position()
#define NEIGHBOR_POSITION_WRITE neighbor.position()



#define NEW_PARTICLE(property_arg) particle.template new_property<property_arg>()




#endif //OPENFPM_PDATA_ALIAS_HPP
