# ===========================================================================
#       
# ===========================================================================
#
# SYNOPSIS
#
#   AX_LIB_HILBERT()
#
# DESCRIPTION
#
#   This macro provides tests of the availability of libHilbert library.
#
#
#   The macro adds a --with-libhilbert option accepting one of three values:
#
#     no   - do not check for the libhilbert library.
#     yes  - do check for libhilbert library in standard locations.
#     path - complete path to the libhilbert library.
#
#   If libhilbert is successfully found, this macro calls
#
#     AC_SUBST(LIBHILBERT_INCLUDE)
#     AC_SUBST(LIBHILBERT_LIB)
#     AC_DEFINE(HAVE_LIBHILBERT)
#
#   and sets with_libhilbert="yes"
#
#   If libhilbert is disabled or not found, this macros sets with_libhilbert="no"
#
#   Your configuration script can test $with_libhilbert to take any further
#   actions. LIBHILBERT_{INCLUDE,LIB} may be used when building with C or C++.
#
#   To use the macro, one would code one of the following in "configure.ac"
#   before AC_OUTPUT:
#
#     1) dnl Check for libhilbert support
#        AX_LIB_HILBERT()
#
#   One could test $with_libhilbert for the outcome or display it as follows
#
#     echo "libhilbert support:  $with_libhilbert"
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

AC_DEFUN([AX_LIB_HILBERT], [
        AC_MSG_CHECKING(for libhilbert library)
        AC_REQUIRE([AC_PROG_CC])
        #
        # User hints...
        #
        AC_ARG_VAR([LIBHILBERT], [Libhilbert library location])
        AC_ARG_WITH([libhilbert],
                [AC_HELP_STRING([--with-libhilbert],
                [user defined path to LIBHILBERT library])],
                [
                        if test -n "$LIBHILBERT" ; then
                                AC_MSG_RESULT(yes)
                                with_libhilbert=$LIBHILBERT
                        elif test "$withval" != no ; then
                                AC_MSG_RESULT(yes)
                                with_libhilbert=$withval
                        else
                                AC_MSG_RESULT(no)
                        fi
                ],
                [
                        if test -n "$PETSC" ; then
                                with_libhilbert=$PETSC
                                AC_MSG_RESULT(yes)
                        else
                                with_petsc=/usr
                                if test ! -f "$with_libhilbert/include/hilbertKey.h" ; then
                                        with_libhilbert=/usr/local
                                        if test ! -f "$with_libhilbert/include/hilbertKey.h" ; then
                                                with_libhilbert=""
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
        # locate LIBHILBERT library
        #

                if test -n "$with_libhilbert" ; then
                        old_CC=$CC
                        old_CFLAGS=$CFLAGS
                        old_LDFLAGS=$LDFLAGS
                        CFLAGS="-I$with_libhilbert/include "
                        LDFLAGS="-L$with_libhilbert/lib "
			CC=$CXX

                        AC_LANG_SAVE
                        AC_LANG_C

                        AC_CHECK_HEADER([hilbertKey.h],libhilbert_h=yes,## Copy LIB and include in the target directory
AC_MSG_WARN([could not find header file hilbertKey.h]))
                        AC_CHECK_LIB([libhilbert],[getIntCoordFromHKey],libhilbert_lib=yes,AC_MSG_WARN([could not find libhilbert]))
                        
                        AC_LANG_RESTORE

                        CFLAGS=$old_CFLAGS
                        LDFLAGS=$old_LDFLAGS
                        CC=$old_CC

                        AC_MSG_CHECKING(LIBHILBERT in $with_libhilbert)
                        if test x"$libhilbert_lib" = x"yes" -a x"$libhilbert_h" = x"yes" ; then
                                AC_SUBST(LIBHILBERT_INCLUDE, [-I$with_libhilbert/include])
                                AC_SUBST(LIBHILBERT_LIB, ["-L$with_libhilbert/lib -llibhilbert"])
                                AC_MSG_RESULT(ok)
 				AC_DEFINE(HAVE_LIBHILBERT,1,[Define if you have LIBHILBERT library])
                        else
                                AC_MSG_RESULT(failed)
                        fi
                fi
                #
                #
                #
                if test x = x"$LIBHILBERT_LIB" ; then
                        ifelse([$2],,[],[$2])
                        :
                else
                        ifelse([$1],,[],[$1])
                        :
                fi
        ])dnl AX_LIB_HILBERT
        
