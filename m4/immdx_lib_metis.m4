AC_DEFUN([IMMDX_LIB_METIS], [
        AC_MSG_CHECKING(for METIS library)
        AC_REQUIRE([AC_PROG_CC])
        #
        # User hints...
        #
        AC_ARG_VAR([METIS], [METIS library location])
        AC_ARG_WITH([metis],
                [AC_HELP_STRING([--with-metis],
                [user defined path to METIS library])],
                [
                        if test -n "$METIS" ; then
                                AC_MSG_RESULT(yes)
                                with_metis=$METIS
                        elif test "$withval" != no ; then
                                AC_MSG_RESULT(yes)
                                with_metis=$withval
                        else
                                AC_MSG_RESULT(no)
                        fi
                ],
                [
                        if test -n "$METIS" ; then
                                with_metis=$METIS
                                AC_MSG_RESULT(yes)
                        else
                                with_metis=/usr
                                if test ! -f "$with_metis/include/metis.h" ; then
                                        with_metis=/usr/local
                                        if test ! -f "$with_metis/include/metis.h" ; then
                                                with_metis=""
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
        # locate METIS library
        #
                if test -n "$with_metis" ; then
                        old_CFLAGS=$CFLAGS
                        old_LDFLAGS=$LDFLAGS
                        CFLAGS="-I$with_metis/include"
                        LDFLAGS="-L$with_metis/lib"

                        AC_LANG_SAVE
                        AC_LANG_C

                        AC_CHECK_LIB(metis, METIS_PartMeshDual,
                                [metis_lib=yes], [metis_lib=yes], [-lm])
                        AC_CHECK_HEADER(metis.h, [metis_h=yes],
                                [metis_h=no], [/* check */])

                        AC_LANG_RESTORE

                        CFLAGS=$old_CFLAGS
                        LDFLAGS=$old_LDFLAGS

                        AC_MSG_CHECKING(METIS in $with_metis)
                        if test "$metis_lib" = "yes" -a "$metis_h" = "yes" ; then
                                AC_SUBST(METIS_INCLUDE, [-I$with_metis/include])
                                AC_SUBST(METIS_LIB, [-L$with_metis/lib])
                                AC_MSG_RESULT(ok)
                        else
                                AC_MSG_RESULT(failed)
                        fi
                fi
                #
                #
                #
                if test x = x"$METIS_LIB" ; then
                        ifelse([$2],,[AC_MSG_ERROR(Failed to find valid METIS library)],[$2])
                        :
                else
                        ifelse([$1],,[AC_DEFINE(HAVE_METIS,1,[Define if you have METIS library])],[$1])
                        :
                fi
        ])dnl IMMDX_LIB_METIS
