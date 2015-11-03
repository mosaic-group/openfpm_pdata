AC_DEFUN([IMMDX_LIB_PARMETIS], [
        AC_MSG_CHECKING(for PARMETIS library)
        AC_REQUIRE([AC_PROG_CC])
        #
        # User hints...
        #
        AC_ARG_VAR([PARMETIS], [PARMETIS library location])
        AC_ARG_WITH([parmetis],
                [AC_HELP_STRING([--with-parmetis],
                [user defined path to PARMETIS library])],
                [
                        if test -n "$PARMETIS" ; then
                                AC_MSG_RESULT(yes)
                                with_parmetis=$PARMETIS
                        elif test "$withval" != no ; then
                                AC_MSG_RESULT(yes)
                                with_parmetis=$withval
                        else
                                AC_MSG_RESULT(no)
                        fi
                ],
                [
                        if test -n "$PARMETIS" ; then
                                with_parmetis=$PARMETIS
                                AC_MSG_RESULT(yes)
                        else
                                with_parmetis=/usr
                                if test ! -f "$with_parmetis/include/parmetis.h" ; then
                                        with_parmetis=/usr/local
                                        if test ! -f "$with_parmetis/include/parmetis.h" ; then
                                                with_parmetis=""
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
        # locate PARMETIS library
        #
                if test -n "$with_parmetis" ; then
                        old_CC=$CC
                        old_CFLAGS=$CFLAGS
                        old_LDFLAGS=$LDFLAGS
                        CFLAGS="-I$with_parmetis/include -I$with_metis/include"
                        LDFLAGS="-L$with_parmetis/lib -L$with_metis/lib"
			CC=$CXX

                        AC_LANG_SAVE
                        AC_LANG_C

                        AC_CHECK_LIB(parmetis, ParMETIS_V3_PartKway,
                                [parmetis_lib=yes], [parmetis_lib=yes], [-lm])
                        AC_CHECK_HEADER(parmetis.h, [parmetis_h=yes],
                                [parmetis_h=no], [/* check */])

                        AC_LANG_RESTORE

                        CFLAGS=$old_CFLAGS
                        LDFLAGS=$old_LDFLAGS
                        CC=$old_CC

                        AC_MSG_CHECKING(PARMETIS in $with_parmetis)
                        if test "$parmetis_lib" = "yes" -a "$parmetis_h" = "yes" ; then
                                AC_SUBST(PARMETIS_INCLUDE, [-I$with_parmetis/include])
                                AC_SUBST(PARMETIS_LIB, [-L$with_parmetis/lib])
                                AC_MSG_RESULT(ok)
                        else
                                AC_MSG_RESULT(failed)
                        fi
                fi
                #
                #
                #
                if test x = x"$PARMETIS_LIB" ; then
                        ifelse([$2],,[AC_MSG_ERROR(Failed to find valid PARMETIS library)],[$2])
                        :
                else
                        ifelse([$1],,[AC_DEFINE(HAVE_PARMETIS,1,[Define if you have PARMETIS library])],[$1])
                        :
                fi
        ])dnl IMMDX_LIB_PARMETIS
