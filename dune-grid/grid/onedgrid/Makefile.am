# $Id: Makefile.am 3085 2006-09-08 09:24:23Z mblatt $

noinst_LTLIBRARIES = libonedgrid.la

libonedgrid_la_SOURCES = onedgrid.cc
libonedgrid_la_LIBADD = $(DUNE_COMMON_LIBS)

onedgriddir = $(includedir)/dune/grid/onedgrid/
onedgrid_HEADERS = \
   onedgridgeometry.hh onedgridentity.hh onedgridhieriterator.hh \
   onedgridleveliterator.hh onedintersectionit.hh onedgridentitypointer.hh \
   onedgridindexsets.hh onedgridleafiterator.hh

include $(top_srcdir)/am/global-rules
