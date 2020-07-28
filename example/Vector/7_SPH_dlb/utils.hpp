#ifndef OPENFPM_PDATA_UTILS_HPP
#define OPENFPM_PDATA_UTILS_HPP

bool amIMaster(Vcluster<>& v_cl) { return v_cl.getProcessUnitID() == 0; }

void printMe(Vcluster<>& v_cl) {
  std::cout << "VCL #" << v_cl.getProcessUnitID() << " ";
}

#endif  // OPENFPM_PDATA_UTILS_HPP
