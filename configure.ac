AC_INIT([BeatBox], [1.7], [V.N.Biktashev@exeter.ac.uk])
AC_CONFIG_AUX_DIR([.])
AM_INIT_AUTOMAKE([-Wall gnu])
myhost=`uname -n`
case $myhost in
   hector*) 
     LIBS=" -lm -lX11 -lxcb -lxcb-xlib -ldl -lXau";;
   *) 
   ;;
   esac
myos=`uname`
case $myos in
   Darwin*)
     LIBS="-lX11 -L/opt/X11/lib"
     CFLAGS="$CFLAGS -I/opt/X11/include"
   ;;
   *)
   ;;
   esac
# echo "#define SUBVERSION \""`svnversion`"\"" > subversion.h  # this does not achieve the purpose as applies outside svn
# AC_DEFINE_UNQUOTED(SUBVERSION, "`svnversion`", [subversion number]) # this does not achieve the purpose as applies outside svn
AC_SEARCH_LIBS(XOpenDisplay, [X11], [AC_SUBST(X11_AVAIL,yes)], [AC_MSG_WARN([X11 libraries not found, Beatbox will be built without graphical displays.  If you need graphical displays you will need to configure Beatbox with X11, please contact the Beatbox developers for further information.])
 AC_SUBST(X11_AVAIL,no)]
, [])
AC_SEARCH_LIBS(glXQueryExtension, [GL], [AC_SUBST(GL_AVAIL,yes)], [AC_MSG_WARN([OpenGL libraries not found, Beatbox will be built without GL graphical extention. If you need graphical displays you will need to configure Beatbox with X11, please contact the Beatbox developers for further information.])
 AC_SUBST(GL_AVAIL,no)]
, [])
AC_SEARCH_LIBS(pow, [m], [], [])
AC_PROG_CC(mpicc gcc cl cc)
AM_PROG_CC_C_O
AC_CONFIG_HEADERS([config.h])
AC_TRY_LINK([#include "mpi.h"],[MPI_Init(0, 0);], [mpi_link="yes"], [mpi_link="no"])
AC_CONFIG_FILES([
	Makefile
	src/Makefile
])

##########################
# As instructed in 
# http://www.mhoenicka.de/system-cgi/blog/index.php?itemid=999

# see whether subversion is installed
AC_PATH_PROG(svnversioncommand, svnversion)

# use svnversion to record the current repository revision only if
# subversion is installed and we are in a working copy
if test "X$svnversioncommand" = "X" || test `$svnversioncommand -n '.'` = "exported"; then
  mysvnversion="cat svn_dist_version"
else
  mysvnversion="svnversion -n '.'"
fi
AM_CONDITIONAL(HAVESVNWC, test "X$svnversioncommand" != "X" && test `$svnversioncommand -n '.'` != "exported")
AC_SUBST(mysvnversion)
##########################

if test "$X11_AVAIL" = "no"; then
  echo "***************"
  echo "Reminder: X11 libraries not found so the graphical part of Beatbox is not being built.  If you want graphical display then please contact the Beatbox developers for help building with X11 libraries."
  echo "***************"
  CFLAGS="$CFLAGS -D NOX11"
else
  CFLAGS="$CFLAGS -D X11"
fi

if test "$GL_AVAIL" = "no"; then
  echo "***************"
  echo "Reminder: OpenGL libraries not found, Beatbox will be built without GL graphical extention. If you need graphical displays you will need to configure Beatbox with X11, please contact the Beatbox developers for further information."
  echo "***************"
  CFLAGS="$CFLAGS -D NOGL"
else
  CFLAGS="$CFLAGS -D GL"
fi

AC_CONFIG_SUBDIRS([src/gsl-1.16.extract])
AC_OUTPUT
