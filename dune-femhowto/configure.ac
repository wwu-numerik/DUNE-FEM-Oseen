# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.50)
DUNE_AC_INIT
AM_INIT_AUTOMAKE([foreign 1.9 tar-pax dist-bzip2])


# check if we're in the right directory
AC_CONFIG_SRCDIR([src_grid/dune_hello_world.cc])
# which header-file is created
AM_CONFIG_HEADER(config.h)

# we need no more than the standard DUNE-stuff
DUNE_CHECK_ALL_M

# if istl not found, its disabled 
#DUNE_CHECK_OPTIONAL([dune-istl],[istl/allocator.hh])

# check dennis package
DUNE_PATH_PARDG

DUNE_SUMMARY_ALL

if test x$HAVE_DUNE_COMMON != x1 ; then
  AC_MSG_ERROR([Can't work without the dune-common module. Maybe you have to supply your dune-common directory as --with-dunecommon=dir])
fi

if test x$HAVE_DUNE_GRID != x1 ; then
  AC_MSG_ERROR([Can't work without the dune-grid module. Maybe you have to supply your dune-common directory as --with-dunegrid=dir])
fi

if test x$HAVE_DUNE_FEM != x1 ; then
  AC_MSG_ERROR([Can't work without the dune-fem module. Maybe you have to supply your dune-common directory as --with-dunegrid=dir])
fi

# preset variable to path such that #include <dune/...> works
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS) -I$(top_srcdir)')
AC_SUBST(AM_LDFLAGS, '$(DUNE_LDFLAGS) $(DUNE_LIBS)')

# implicitly set the Dune-flags everywhere
CPPFLAGS="$DUNE_CPPFLAGS"
LDFLAGS="$DUNE_LDFLAGS"
LIBS="$DUNE_LIBS $BLAS_LIBS"

# generate Makefiles
AC_CONFIG_FILES([Makefile
     src_grid/Makefile
     src_poisson/Makefile
     src_finitevolume/Makefile	
     src_localdg/Makefile
     sources/Makefile
     sources/poisson/Makefile
     sources/localdg/Makefile
     sources/localdg/advectiondiffusion/Makefile
     ])
AC_OUTPUT
