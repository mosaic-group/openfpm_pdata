# ===========================================================================
#        http://www.gnu.org/software/autoconf-archive/ax_lib_petsc.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_PETSC()
#
# DESCRIPTION
#
#   This macro provides tests of the availability of PETSC library.
#
#
#   The macro adds a --with-petsc option accepting one of three values:
#
#     no   - do not check for the PETSC library.
#     yes  - do check for PETSC library in standard locations.
#     path - complete path to the PETSC library.
#
#   If PETSC is successfully found, this macro calls
#
#     AC_SUBST(PETSC_INCLUDE)
#     AC_SUBST(PETSC_LIB)
#     AC_DEFINE(HAVE_PETSC)
#
#   and sets with_petsc="yes"
#
#   If PETSC is disabled or not found, this macros sets with_petsc="no"
#
#   Your configuration script can test $with_petsc to take any further
#   actions. PETSC_{INCLUDE,LIB} may be used when building with C or C++.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for PETSC support
#        AX_LIB_PETSC()
#
#   One could test $with_petsc for the outcome or display it as follows
#
#     echo "PETSC support:  $with_petsc"
#
#   You could also for example, override the default CC in "configure.ac"
#
# LICENSE
#
#   Copyright (c) 2009 Timothy Brown <tbrown@freeshell.org>
#   Copyright (c) 2010 Rhys Ulerich <rhys.ulerich@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 12

AC_DEFUN([AX_LIB_PETSC], [
        AC_MSG_CHECKING(for PETSC library)
        AC_REQUIRE([AC_PROG_CC])
        #
        # User hints...
        #
        AC_ARG_VAR([PETSC], [PETSC library location])
        AC_ARG_WITH([petsc],
                [AC_HELP_STRING([--with-petsc],
                [user defined path to PETSC library])],
                [
                        if test -n "$PETSC" ; then
                                AC_MSG_RESULT(yes)
                                with_petsc=$PETSC
                        elif test "$withval" != no ; then
                                AC_MSG_RESULT(yes)
                                with_petsc=$withval
                        else
                                AC_MSG_RESULT(no)
                        fi
                ],
                [
                        if test -n "$PETSC" ; then
                                with_petsc=$PETSC
                                AC_MSG_RESULT(yes)
                        else
                                with_petsc=/usr
                                if test ! -f "$with_petsc/include/petsc.h" ; then
                                        with_petsc=/usr/local
                                        if test ! -f "$with_petsc/include/petsc.h" ; then
                                                with_petsc=""
                                                AC_MSG_RESULT(failed)
                                        else
                                                AC_MSG_RESULT(yes)
                                        fi
                                else
                                        AC_MSG_RESULT(yes)
                                fi
                        fi
                ])
        #
        # locate PETSC library
        #

                if test -n "$with_petsc" ; then
                        old_CC=$CC
                        old_CFLAGS=$CFLAGS
                        old_LDFLAGS=$LDFLAGS
                        CFLAGS="-I$with_petsc/include $HDF5_INCLUDE $METIS_INCLUDE "
                        LDFLAGS="-L$with_petsc/lib $HDF5_LDFLAGS  $HDF5_LIBS $METIS_LIB -lmetis "
			CC=$CXX

                        AC_LANG_SAVE
                        AC_LANG_C

                        AC_CHECK_HEADER([petsc.h],petsc_h=yes,## Copy LIB and include in the target directory
AC_MSG_WARN([could not find header file petsc.h]))
                        AC_CHECK_LIB([petsc],[PetscTrMalloc],petsc_lib=yes,AC_MSG_WARN([could not find libpetsc]))
                                
                        AC_LANG_RESTORE

                        CFLAGS=$old_CFLAGS
                        LDFLAGS=$old_LDFLAGS
                        CC=$old_CC

                        AC_MSG_CHECKING(PETSC in $with_petsc)
                        if test x"$petsc_lib" = x"yes" -a x"$petsc_h" = x"yes" ; then
                                AC_SUBST(PETSC_INCLUDE, [-I$with_petsc/include])
                                AC_SUBST(PETSC_LIB, ["-L$with_petsc/lib -lpetsc"])
                                AC_MSG_RESULT(ok)
 				AC_DEFINE(HAVE_PETSC,1,[Define if you have PETSC library])
                        else
                                AC_MSG_RESULT(failed)
                        fi
                fi
                #
                #
                #
                if test x = x"$PETSC_LIB" ; then
                        ifelse([$2],,[],[$2])
                        :
                else
                        ifelse([$1],,[],[$1])
                        :
                fi
        ])dnl AX_LIB_PETSC
