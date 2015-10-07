#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

## Take all the options with the exception of --enable-install-req

AC_PREREQ(2.59)
AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_CANONICAL_SYSTEM
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_SUBDIRS([openfpm_data openfpm_devices openfpm_vcluster openfpm_io openfpm_numerics])

#### Adding --with-pdata option and openfpm_pdata to prefix folder

if test "$prefix" = "NONE"; then
    prefix=/usr/local
fi
base=$prefix
prefix="$prefix/openfpm_pdata"
echo "Installation dir is: $prefix"

ac_configure_args="$ac_configure_args --with-pdata=../../src"

########################

AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([src/config/config.h])
m4_ifdef([ACX_PTHREAD],,[m4_include([m4/acx_pthread.m4])])
m4_ifdef([AX_BOOST],,[m4_include([m4/ax_boost.m4])])
m4_ifdef([ACX_MPI],,[m4_include([m4/acx_mpi.m4])])
m4_ifdef([AX_OPENMP],,[m4_include([m4/ax_openmp.m4])])
m4_ifdef([AX_CUDA],,[m4_include([m4/ax_cuda.m4])])
m4_ifdef([IMMDX_LIB_METIS],,[m4_include([m4/immdx_lib_metis.m4])])

CXXFLAGS+=" --std=c++11 "
NVCCFLAGS=" "
INCLUDES_PATH=" "

# Create a file with the install base folder
echo "$base" > install_dir

# Checks for programs.
AC_PROG_CXX

# Checks g++ flags

AC_CANONICAL_HOST

# Check that the compiler support mpi

AC_LANG_PUSH([C++])
AC_CHECK_HEADER(mpi.h,[],[echo "mpi.h not found"
                          exit 200])
AC_LANG_POP([C++])

# Check target architetture

#AX_GCC_ARCHFLAG([], [CXXFLAGS="$CXXFLAGS $ax_cv_gcc_archflag"], [])

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
	CXXFLAGS="$CXXFLAGS -Wall -O3 -g3 -funroll-loops "
	NVCCFLAGS+="$NVCCFLAGS -O3 "
fi

## Check for memcheck

AC_MSG_CHECKING(whether to build with memcheck capabilities)
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

## Check for Metis

IMMDX_LIB_METIS([],[echo "Cannot detect metis, use the --with-metis option if it is not installed in the default location"
                    exit 201])

####### include OpenFPM_devices include path

INCLUDES_PATH+="-I. -Isrc/config/ -I../openfpm_io/src -I../openfpm_data/src -I../openfpm_devices/src -I../openfpm_vcluster/src/"

###### Check for se-class1

AC_MSG_CHECKING(whether to build with security enhancement class1)
se_class1=no
AC_ARG_ENABLE(se-class1,
        AC_HELP_STRING(
            [--enable-se-class1],
            [enable security enhancement class1]
        ),
        se_class1="$enableval"
)



AC_MSG_RESULT($se_class1)
if test x"$se_class1" = x"yes"; then
        AC_DEFINE([SE_CLASS1],[],[Security enhancement class 1])
fi


###### Check for se-class 2

AC_MSG_CHECKING(whether to build with security enhancement class 2)
se_class2=no
AC_ARG_ENABLE(se-class2,
        AC_HELP_STRING(
            [--enable-se-class2],
            [enable security enhancement class 2]
        ),
        se_class2="$enableval"
)


AC_MSG_RESULT($se_class2)
if test x"$se_class2" = x"yes"; then
        AC_DEFINE([SE_CLASS2],[],[Security enhancement class 2])
fi

###### Check for se-class 3

AC_MSG_CHECKING(whether to build with security enhancement class 3)
se_class3=no
AC_ARG_ENABLE(se-class3,
        AC_HELP_STRING(
            [--enable-se-class3],
            [enable security enhancement class 3]
        ),
        se_class3="$enableval"
)



AC_MSG_RESULT($se_class3)
if test x"$se_class3" = x"yes"; then
        AC_DEFINE([SE_CLASS3],[],[Security enhancement class 3])
fi

###### Check for action on error

action_on_e=continue
AC_ARG_WITH([action-on-error],
             AS_HELP_STRING([--with-action-on-error=stop,throw,continue],
             [specify the action to do in case of error]),
             [action_on_e="$withval"],
             [action_on_e=continue])

if test x"$action_on_e" = x"stop"; then
           AC_DEFINE([STOP_ON_ERROR],[],[If an error occur stop the program])
fi
if test x"$action_on_e" = x"throw"; then
           AC_DEFINE([THROW_ON_ERROR],[],[when an error accur continue but avoid unsafe operation])
fi

##### CHECK FOR BOOST ##############

AX_BOOST([1.52],[],[echo "boost not found"
                    exit 202])

####### Checking for GPU support

AX_CUDA

## detect for NVCC

if test x"$NVCC_EXIST" = x"yes"; then
  AC_MSG_CHECKING(whether to build with GPU support)
  gpu_support=yes
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
  else
	CUDA_LIBS=""
	CUDA_CFLAGS=""
  fi
else
  gpu_support=no
fi

if test x$gpu_support = x"no"; then
     CUDA_LIBS=""
     CUDA_CFLAGS=""
fi

# Set this conditional if cuda is wanted

AM_CONDITIONAL(BUILDCUDA, test ! x$NVCC = x"no")

###########################

# Define that there is MPI

AC_DEFINE([HAVE_MPI],[],[MPI Enabled])

AC_SUBST(NVCCFLAGS)
AC_SUBST(INCLUDES_PATH)

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile vtk/Makefile ])
AC_OUTPUT
echo ""
echo   "***********************************"
echo   "*                                 *"
if [ test x"$profiler" = x"yes" ]; then
  echo "*    profiler: yes                *"
else
  echo "*    profiler: no                 *"
fi
if [ test x"$memcheck" = x"yes" ]; then
  echo "*    memcheck: yes                *"
else
  echo "*    memcheck: no                 *"
fi
if [ test x"$debuger" = x"yes" ]; then
  echo "*    debug: yes                   *"
else
  echo "*    debug: no                    *"
fi
if [ test x"$se_class1" = x"yes" ]; then
  echo "*    se-class1: yes               *"
else
  echo "*    se-class1: no                *"
fi
if [ test x"$se_class2" = x"yes" ]; then
  echo "*    se-class2: yes               *"
else
  echo "*    se-class2: no                *"
fi
if [ test x"$se_class3" = x"yes" ]; then
  echo "*    se-class3: yes               *"
else
  echo "*    se-class3: no                *"
fi
if [ test x"$gpu_support" = x"no" ]; then
  echo "*    gpu: no                      *"
else
  echo "*    gpu: yes                     *"
fi
echo   "*                                 *"
echo   "***********************************"

