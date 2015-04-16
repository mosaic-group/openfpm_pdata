#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_PREREQ(2.59)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CANONICAL_SYSTEM
AC_CONFIG_SRCDIR([src/main.cpp])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([src/config/config.h])
m4_ifdef([MYSQL_FOUND],,[m4_include([m4/ax_lib_mysql.m4])])
m4_ifdef([AX_CHECK_COMPILER_FLAGS],,[m4_include([m4/ax_check_compiler_flags.m4])])
m4_ifdef([ACX_PTHREAD],,[m4_include([m4/acx_pthread.m4])])
m4_ifdef([AX_CHECK_CL],,[m4_include([m4/ax_opencl.m4])])
m4_ifdef([AX_BOOST_BASE],,[m4_include([m4/ax_boost_base.m4])])
m4_ifdef([AX_BOOST_PROGRAM_OPTIONS],,[m4_include([m4/ax_boost_program_options.m4])])
m4_ifdef([AX_BOOST_THREAD],,[m4_include([m4/ax_boost_thread.m4])])
m4_ifdef([ACX_MPI],,[m4_include([m4/acx_mpi.m4])])
m4_ifdef([AX_OPENMP],,[m4_include([m4/ax_openmp.m4])])
m4_ifdef([AX_GCC_X86_CPUID],,[m4_include([m4/ax_gcc_x86_cpuid.m4])])
m4_ifdef([AX_GCC_ARCHFLAG],,[m4_include([m4/ax_gcc_archflag.m4])])
m4_ifdef([AX_CUDA],,[m4_include([m4/ax_cuda.m4])])

CXXFLAGS+=" --std=c++11 "
NVCCFLAGS=" "
INCLUDES_PATH=" "

# Checks for programs.
AC_PROG_CXX

# Checks g++ flags

AC_CANONICAL_HOST

# Check target architetture

AX_GCC_ARCHFLAG([], [CXXFLAGS="$CXXFLAGS $ax_cv_gcc_archflag"], [])

###### Check for debug compilation

AC_MSG_CHECKING(whether to build with debug information)
debuger=no
AC_ARG_ENABLE(debug,
	AC_HELP_STRING(
	    [--enable-debug],
	    [enable debug data generation (def=no)]
	),
  	debuger="$enableval"
)



AC_MSG_RESULT($debuger)
if test x"$debuger" = x"yes"; then
	AC_DEFINE([DEBUG_MODE],[],[Debug])
	AC_DEFINE([DEBUG],[],[Debug])
	CXXFLAGS="$CXXFLAGS -g3  -Wall -O0 "
	NVCCFLAGS+="$NVCCFLAGS -g -O0 "
else
	CXXFLAGS="$CXXFLAGS -Wall -O3 -g3 "
	NVCCFLAGS+="$NVCCFLAGS -O3 "
fi

## Check for memcheck

AC_MSG_CHECKING(whether to build with memcheck capabilities)
debuger=no
AC_ARG_ENABLE(memcheck,
        AC_HELP_STRING(
            [--enable-memcheck],
            [enable memory check (def=no)]
        ),
        memcheck="$enableval"
)

AC_MSG_RESULT($memcheck)

if test x"$memcheck" = x"yes"; then
        AC_DEFINE([MEMLEAK_CHECK],[],[Memory check, corruption and leak])
fi

#########

####### include OpenFPM_devices include path

INCLUDES_PATH+="-I. -Isrc/config/ -I../../metis_install/include -I../../OpenFPM_IO/src -I../../OpenFPM_data/src -I../../OpenFPM_devices/src -I../../OpenFPM_vcluster/src/"

####### Checking for GPU support

AX_CUDA

AC_MSG_CHECKING(whether to build with GPU support)
gpu_support=no
AC_ARG_ENABLE(gpu,
	AC_HELP_STRING(
	    [--enable-gpu],
	    [enable gpu support]
	),
  	gpu_support="$enableval"
)



AC_MSG_RESULT($gpu_support)
if test x"$gpu_support" = x"yes"; then
	AC_DEFINE([GPU],[],[GPU support])
fi


# Set this conditional if cuda is wanted

AM_CONDITIONAL(BUILDCUDA, test ! x$NVCC = x"no")

###########################

# Define that there is MPI

AC_DEFINE([HAVE_MPI],[],[MPI Enabled])

#####

no_avx=no
no_sse42=no
no_sse41=no
no_sse3=no
no_sse2=no
no_sse=no
no_mmx=no
AX_CHECK_COMPILER_FLAGS([-msse4.2],[CXXFLAGS="$CXXFLAGS -mavx"],[no_avx=yes])
AX_CHECK_COMPILER_FLAGS([-msse4.2],[CXXFLAGS="$CXXFLAGS -msse4.2"],[no_sse42=yes])
AX_CHECK_COMPILER_FLAGS([-msse4.1],[CXXFLAGS="$CXXFLAGS -msse4.1"],[no_sse41=yes])
AX_CHECK_COMPILER_FLAGS([-msse3],[CXXFLAGS="$CXXFLAGS -msse3"],[no_sse3=yes])
AX_CHECK_COMPILER_FLAGS([-msse2],[CXXFLAGS="$CXXFLAGS -msse2"],[no_sse2=yes])
AX_CHECK_COMPILER_FLAGS([-msse],[CXXFLAGS="$CXXFLAGS -msse"],[no_sse=yes])
AX_CHECK_COMPILER_FLAGS([-mmmx],[CXXFLAGS="$CXXFLAGS -mmmx"],[no_mmx=yes])
AX_CHECK_COMPILER_FLAGS([-Wno-unused-but-set-variable],[CXXFLAGS="$CXXFLAGS -Wno-unused-but-set-variable"],[])

AC_SUBST(NVCCFLAGS)
AC_SUBST(INCLUDES_PATH)

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
echo ""
echo   "***********************************"
echo   "*                                 *"
arch_str="${ax_cv_gcc_archflag#-march=#-mtune=}"
arch_str="${arch_str#-mtune=}"
n_arch_str=${#arch_str}
for (( X=0; X<23-n_arch_str; X++ ))
do
   arch_str="$arch_str "
done
echo "*    arch: $arch_str*"
if [ test x"$no_sse42" = x"no" ]; then
    echo "*    sse4.2: yes                  *"
else
    echo "*    sse4.2: no                   *"
fi
if [ test x"$no_sse41" = x"no" ]; then
    echo "*    sse4.1: yes                  *"
else
    echo "*    sse4.1: no                   *"
fi
if [ test x"$no_sse3" = x"no" ]; then
    echo "*    sse3: yes                    *"
else
    echo "*    sse3: no                     *"
fi
if [ test x"$no_sse2" = x"no" ]; then
    echo "*    sse2: yes                    *"
else
    echo "*    sse2: no                     *"
fi
if [ test x"$no_sse" = x"no" ]; then
    echo "*    sse: yes                     *"
else
    echo "*    sse: no                      *"
fi
if [ test x"$no_mmx" = x"no" ]; then
  echo "*    mmx: yes                     *"
else
  echo "*    mmx: no                      *"
fi
if [ test x"$profiler" = x"yes" ]; then
  echo "*    profiler: yes                *"
else
  echo "*    profiler: no                 *"
fi
if [ test x"$debuger" = x"yes" ]; then
  echo "*    debug: yes                   *"
else
  echo "*    debug: no                    *"
fi
if [ test x"$no_64" = x"no" ]; then
  echo "*    64 bit: yes                  *"
else
  echo "*    64 bit: no                   *"
fi
if [ test x"$gpu_support" = x"no" ]; then
  echo "*    gpu: no                      *"
else
  echo "*    gpu: yes                     *"
fi
echo   "*                                 *"
echo   "***********************************"

