# $Id: xdr.m4 1236 2004-12-08 18:54:46Z thimo $
# searches for XDR Headers which are implicitly included by rpc.h
# some systems don't like it when xdr.h is directly included

AC_DEFUN([DUNE_PATH_XDR],[
  AC_REQUIRE([AC_PROG_CC])
  AC_CHECK_HEADERS(rpc/rpc.h)
])
