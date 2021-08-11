#include "../../src/Vector/vector_dist.hpp"
#include "../../Openfpm_numerics/src/FiniteDifference/Upwind_gradient.hpp"
#include "../../Openfpm_numerics/src/FiniteDifference/Eno_Weno.hpp"
#include "/Library/Developer/CommandLineTools/SDKs/MacOSX11.1.sdk/System/Library/Frameworks/Kernel.framework/Versions/A/Headers/math.h"
#include "../../openfpm_numerics/src/Draw/DrawParticles.hpp"
#include "/Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/System/Library/Frameworks/Kernel.framework/Versions/A/Headers/stddef.h"


#define PHASE_A 0
#define PHASE_B 1

const double dp = 1/64.0;

const size_t type = 0;

int main()
{
	std::cout<<type<<std::endl; // @suppress("Symbol is not resolved")


	return 0;
}
