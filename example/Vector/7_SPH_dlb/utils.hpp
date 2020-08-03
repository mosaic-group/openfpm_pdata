#ifndef OPENFPM_PDATA_UTILS_HPP
#define OPENFPM_PDATA_UTILS_HPP

bool amIMaster(Vcluster<>& v_cl) { return v_cl.getProcessUnitID() == 0; }

void printMe(Vcluster<>& v_cl) {
  std::cout << "VCL #" << v_cl.getProcessUnitID() << " ";
}

template <typename T>
bool isIn(const openfpm::vector<T>& v, const T& x) {
  for (auto i = 0; i < v.size(); i++) {
    if (v.get(i) == x) {
      return true;
    }
  }

  return false;
}

#endif  // OPENFPM_PDATA_UTILS_HPP
