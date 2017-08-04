# ===========================================================================
#         http://www.gnu.org/software/autoconf-archive/ax_lapack.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_SUITESPARSE([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for SuiteSparse library (see http://www.netlib.org/lapack/)
#   On success, it sets the SUITESPARSE_LIBS output variable to hold the 
#   requisite library linkages.
#
#   To link with SUITESPARSE, you should link with:
#
#     $SUITESPARSE_LIBS
#
#   The user may also use --with-suitesparse=<lib> in order to use some specific
#   SuiteSparse library <lib>.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a SUITESPARSE library
#   is found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it
#   is not found. If ACTION-IF-FOUND is not specified, the default action
#   will define HAVE_LAPACK.
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

AU_ALIAS([ACX_SUITESPARSE], [AX_SUITESPARSE])
AC_DEFUN([AX_SUITESPARSE], [
ax_suitesparse_ok=no

AC_ARG_WITH(suitesparse,
        [AS_HELP_STRING([--with-suitesparse=directory], [use SuiteSparse directory])],[
                        SUITESPARSE_LIBS="-L$with_suitesparse/lib"
                        SUITESPARSE_INCLUDE="-I$with_suitesparse/include"
                        ])

#
# Platform specific setup
#
#############################
AC_CANONICAL_HOST
# Check for which host we are on and setup a few things
# specifically based on the host
case $host_os in
  darwin* )
        RT_LIB=" $BLAS_LIBS"
        ;;
  linux*)
        RT_LIB="-lrt $BLAS_LIBS"
        ;;
    *)
        RT_LIB="-lrt $BLAS_LIBS"
        ;;
esac

# First, check SUITESPARSE_LIBS environment variable
if test "x$SUITESPARSE_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$SUITESPARSE_LIBS -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig -lm  $RT_LIB"
        AC_MSG_CHECKING([for umf_l_malloc])
        AC_TRY_LINK_FUNC(umf_l_malloc, [ax_suitesparse_ok=yes
                                        SUITESPARSE_LIBS="$SUITESPARSE_LIBS -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig"
                                       ], [SUITRSPARSE_LIBS=""])
        AC_MSG_RESULT($ax_suitesparse_ok)
        LIBS="$save_LIBS"
        if test $ax_suitesparse_ok = no; then
                SUITESPARSE_LIBS=""
        fi
	old_CFLAGS="$CFLAGS"
        CFLAGS=$SUITESPARSE_INCLUDE
	AC_CHECK_HEADER(umfpack.h,[],[SUITESPARSE_INCLUDE=""
                                     ax_suitesparse_ok=no])
                                     
        CFLAGS="$old_CFLAGS"
else
        AC_CHECK_LIB(umfpack,umf_l_alloc,[SUITESPARSE_LIBS="$SUITESPARSE_LIBS -lumfpack -lamd -lbtf -lcamd -lccolamd -lcholmod -lcolamd -lcxsparse -lklu -ldl -lrbio -lspqr -lsuitesparseconfig"],[
                                                                  SUITESPARSE_LIBS=""
                                                                  ax_suitesparse_ok=no
                                                                  ])
        old_CFLAGS="$CFLAGS"
        AC_CHECK_HEADER(umfpack.h,[],[SUITESPARSE_INCLUDE=""
                                      ax_suitesparse_ok=no])
        
                                      
        CFLAGS="$old_CFLAGS"
fi

AC_SUBST(SUITESPARSE_LIBS)
AC_SUBST(SUITESPARSE_INCLUDE)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_suitesparse_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_SUITESPARSE,1,[Define if you have SUITESPARSE library.]),[$1])
        :
else
        ax_suitesparse_ok=no
        $2
fi
])dnl AX_SUITESPARSE

