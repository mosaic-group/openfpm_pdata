//
// Created by landfried on 29.09.22.
//

#ifndef OPENFPM_PDATA_TP_INSTANCE_HPP
#define OPENFPM_PDATA_TP_INSTANCE_HPP

#include "TP_Algorithm.hpp"


double GlobalVariable::dt = 0.05;
double GlobalVariable::t = 0;
double GlobalVariable::t_final = 10.1;

double GlobalVariable::domainSize = 40.0;
int GlobalVariable::meshSize = 256;
double GlobalVariable::meshSpacing = domainSize / meshSize;






#endif //OPENFPM_PDATA_TP_INSTANCE_HPP
