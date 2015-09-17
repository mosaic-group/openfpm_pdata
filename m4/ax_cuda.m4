# -*- mode: autoconf -*-
#
# AX_OPENCL
#
# Check for an OpenCL implementation. If CL is found, _OPENCL is defined and
# the required compiler and linker flags are included in the output variables
# "CL_CFLAGS" and "CL_LIBS", respectively. If no usable CL implementation is
# found, "no_cl" is set to "yes".
#
# If the header "CL/OpenCL.h" is found, "HAVE_CL_OPENCL_H" is defined. If the
# header "OpenCL/OpenCL.h" is found, HAVE_OPENCL_OPENCL_H is defined. These
# preprocessor definitions may not be mutually exclusive.
#
# Based on AX_CHECK_GL, version: 2.4 author: Braden McDaniel
# <braden@endoframe.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#
# As a special exception, the you may copy, distribute and modify the
# configure scripts that are the output of Autoconf when processing
# the Macro. You need not follow the terms of the GNU General Public
# License when using or distributing such scripts.
#
AC_DEFUN([AX_CUDA],
[AC_REQUIRE([AC_CANONICAL_HOST])dnl

# Search nvcc compiler
AC_CHECK_PROG([NVCC_EXIST],[nvcc],["yes"],["no"])
AS_IF([test "x$NVCC_EXIST" = "xno"],[],[
          NVCC=`which nvcc`

          # Set CUDA_CFLAGS to $NVCC, where substring "bin/nvcc"
          # is substituted by "include".
          CUDA_CFLAGS=" ${NVCC%bin//nvcc}"
          CUDA_CFLAGS=" ${CUDA_CFLAGS%bin/nvcc}"
          CUDA_CFLAGS=" -I${CUDA_CFLAGS}include"

          #Set CUDA_CFLAGS to $NVCC, where substring "bin/nvcc"
          #is substituted by "lib".
          CUDA_LIBS="${NVCC%bin//nvcc}"
          CUDA_LIBS="${CUDA_LIBS%bin/nvcc}"
          CUDA_PATH=$CUDA_LIBS
          CUDA_LIBS=" -L${CUDA_LIBS}lib"

          # If $build_cpu contains "_64", append "64" to CUDA_LIBS
          AS_IF([echo $build_cpu | grep -q "_64"],
                [
                 AS_IF([ test -d $CUDA_PATH/lib64 ], [ CUDA_LIBS+="64" ], [])
                 # Be carefull the return code 0 mean true return code 1 mean false
                 AS_IF([ command -v bumblebee >/dev/null ], [ CUDA_LIBS+=" -L/usr/lib64/nvidia-bumblebee/ "  ],
                                                             [
                                                               echo "bumblebee, NVIDIA optimus,  not found"
                                                             ])
                 AS_IF([ test -d /usr/local/cuda/lib64  ], [ CUDA_LIBS+=" -L/usr/local/cuda/lib64 "  ],
                       [
                        AS_IF([ test -d /usr/local/cuda/lib ],[ CUDA_LIBS+=" -L/usr/local/cuda/lib  "  ])
                       ])
                ])

          # Append " -lcuda -lcudart" to CUDA_LIBS
          CUDA_LIBS+=" -lcuda -lcudart"
          
          # Make variables available in Makefile.am
          AC_SUBST([CUDA_CFLAGS])
          AC_SUBST([CUDA_LIBS])
          echo $NVCC
          AC_SUBST([NVCC])
          AC_DEFINE([NVCC],[],[NVCC compiling])
])dnl

])dnl
