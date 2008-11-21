# $Id: Makefile.am,v 1.15 2006/10/16 15:36:33 robertk Exp $

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# do everything in the subdirs
SUBDIRS = doc src_localdg #sources #only compile things that compile without errors

#pkgconfigdir = $(libdir)/pkgconfig
#pkgconfig_DATA = dune-fem.pc

DISTCHECK_CONFIGURE_FLAGS = --with-dune-common="$(DUNE_COMMON_ROOT)" --with-dune-grid="$(DUNE_GRID_ROOT)" CXX="$(CXX)" CC="$(CC)" DOXYGEN="true"

include $(top_srcdir)/am/global-rules
