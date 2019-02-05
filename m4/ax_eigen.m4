#
#   AX_EIGEN([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for EIGEN library (see http://www.netlib.org/lapack/)
#   On success, it sets the EIGEN_INCLUDE output variable to hold the 
#   requisite includes.
#
#   The user may also use --with-eigen=<include> in order to use some specific
#   Eigen library.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a Eigen library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_EIGEN.
#
# LICENSE
#
#   Copyright (c) 2009 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 1

AU_ALIAS([ACX_EIGEN], [AX_EIGEN])
AC_DEFUN([AX_EIGEN], [
ax_eigen_ok=no

AC_ARG_WITH(eigen,
        [AS_HELP_STRING([--with-eigen=directory], [use Eigen directory])],[
                        EIGEN_INCLUDE="-I$with_eigen/include"
                        ])

# First, check EIGEN_INCLUDE environment variable
if test "x$EIGEN_INCLUDE" != x; then
	old_CXXFLAGS="$CXXFLAGS"
        AC_LANG_PUSH([C++])
	CXXFLAGS+=" -I${withval} -DEIGEN"

        # Check for the EIGEN header files
        AC_CHECK_HEADERS([Eigen/Dense Eigen/LU],[ax_eigen_ok=yes],[ax_eigen_ok=no])
        AC_LANG_POP()
        CFLAGS="$old_CFLAGS"
else
        old_CXXFLAGS="$CXXFLAGS"
        CXXFLAGS+=' -I/usr/include/eigen3 -DEIGEN'
        AC_CHECK_HEADERS([Eigen/Dense Eigen/LU],[ax_eigen_ok=yes],[ax_eigen_ok=no])
        CFLAGS="$old_CFLAGS"
fi

AC_SUBST(EIGEN_INCLUDE)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_eigen_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_EIGEN,1,[Define if you have EIGEN library.]),[$1])
        :
else
        ax_eigen_ok=no
        $2
fi
])dnl AX_EIGEN

