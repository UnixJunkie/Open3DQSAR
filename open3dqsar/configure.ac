# configure.ac
# 
# is part of
# 
# Open3DQSAR
# ----------
# 
# An open-source software aimed at high-throughput
# chemometric analysis of molecular interaction fields
# 
# Copyright (C) 2009-2013 Paolo Tosco, Thomas Balle
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
# 
# For further information, please contact:
# 
# Paolo Tosco, PhD
# Dipartimento di Scienza e Tecnologia del Farmaco
# Universita' degli Studi di Torino
# Via Pietro Giuria, 9
# 10125 Torino (Italy)
# Phone:  +39 011 670 7680
# Mobile: +39 348 553 7206
# Fax:    +39 011 670 7687
# E-mail: paolo.tosco@unito.it

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(Open3DQSAR, 3.0, <paolo.tosco@unito.it>)
AC_DEFINE([PACKAGE_NAME_SMALL_CAPS], ["open3dqsar"], [Package name in small caps])
AC_DEFINE([O3Q], [1], [Open3DQSAR])
AM_INIT_AUTOMAKE(Open3DQSAR, 3.0)
AC_CONFIG_SRCDIR([src])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

# Checks for libraries.
cat << eof > src/test.c
int main(int argc, char **argv)
{
	return 0;
}
eof

AC_CHECK_LIB([m], [sqrt], , AC_MSG_FAILURE([libm not found]))
AC_CHECK_LIB([z], [gzopen], , AC_MSG_FAILURE([libz not found]))
AC_CHECK_LIB([z -lminizip], [zipOpen], , AC_MSG_FAILURE([libminizip not found]))

# Checks for library functions.
AC_FUNC_VPRINTF
AC_FUNC_FORK
AC_TYPE_SIGNAL
AC_CHECK_FUNCS([dup2 getcwd gettimeofday memset mkdir mkdtemp mkstemp munmap \
	pow putenv rint rmdir setenv sqrt strcasecmp strchr strncasecmp strstr strtok_r uname])

LIBEDIT=
AC_CHECK_LIB([edit], [readline], [AC_SUBST([LIBEDIT], ["-ledit"])
	AC_DEFINE([HAVE_EDITLINE_FUNCTIONALITY], [1], [Use editline])])
AS_IF([test -z "$LIBEDIT"],
	[AC_CHECK_LIB([edit -lcurses], [readline], [AC_SUBST([LIBEDIT], ["-ledit -lcurses"])
	AC_DEFINE([HAVE_EDITLINE_FUNCTIONALITY], [1], [Use editline])])])
AC_CHECK_LIB([ws2_32], [main], ,
	AC_CHECK_LIB([pthread], [pthread_create], ,
	AC_MSG_FAILURE([libpthread not found])))
AC_CHECK_LIB([ws2_32], [main],
	AC_CHECK_LIB([shlwapi], [PathMatchSpecA], ,
	AC_MSG_FAILURE([libshlwapi not found])))
AC_CHECK_LIB([ws2_32], [main],
	[AS_IF([$CC $CFLAGS src/test.c src/open3dqsar_win32.res -o src/test 2> /dev/null],
	[AC_SUBST([WINDOWS_RES], ["open3dqsar_win32.res"])],
	[AC_SUBST([WINDOWS_RES], ["open3dqsar_win64.res"])])])
rm -f src/test.c src/test

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h stddef.h stdlib.h string.h sys/param.h sys/time.h termios.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_HEADER_TIME

AC_ARG_WITH([atlas],
	[AC_HELP_STRING([--with-atlas],
	[link against ATLAS CBLAS and LAPACK math libraries])],
	[],
	[with_atlas=no])

LIBATLAS=
LIBLAPACK_ATLAS=
LIBLAPACK_F77=
AS_IF([test "x$with_atlas" != xno],
	[AC_CHECK_HEADER([clapack.h],
	AC_DEFINE([HAVE_CLAPACK_H], [1], [Use clapack.h header]),
	[AC_MSG_FAILURE(
	[--with-atlas was given, but test for clapack.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_atlas" != xno],
	[AC_CHECK_HEADER([cblas.h],
	AC_DEFINE([HAVE_CBLAS_H], [1], [Use cblas.h header]),
	[AC_MSG_FAILURE(
	[--with-atlas was given, but test for cblas.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_atlas" != xno],
	[AC_CHECK_LIB([cblas -latlas], [cblas_dgemm],
	[AC_SUBST([LIBATLAS], ["-lcblas -latlas"])
	AC_DEFINE([HAVE_LIBATLAS], [1], [Use ATLAS CBLAS library])])])
AS_IF([test "x$with_atlas" != xno],
	[AC_CHECK_LIB([lapack -lcblas -latlas], [clapack_dgetrf],
	[AC_SUBST([LIBLAPACK_ATLAS], ["-llapack"])
	AC_DEFINE([HAVE_LIBLAPACK_ATLAS], [1], [Use ATLAS LAPACK library])])])
AS_IF([test "x$with_atlas" != xno && test -z "$LIBLAPACK_ATLAS"],
	[AC_CHECK_LIB([lapack -lalapack -lcblas -latlas], [clapack_dgetrf],
	[AC_SUBST([LIBLAPACK_ATLAS], ["-llapack -lalapack"])
	AC_DEFINE([HAVE_LIBLAPACK_ATLAS], [1], [Use ATLAS LAPACK library])])])
AS_IF([test "x$with_atlas" != xno && test -z "$LIBLAPACK_ATLAS"],
	[AC_CHECK_LIB([lapack -lf77blas -lcblas -latlas], [clapack_dgetrf],
	[AC_SUBST([LIBLAPACK_ATLAS], ["-llapack -lf77blas"])
	AC_DEFINE([HAVE_LIBLAPACK_ATLAS], [1], [Use ATLAS LAPACK library])])])
AS_IF([test "x$with_atlas" != xno && test -z "$LIBLAPACK_ATLAS"],
	[AC_CHECK_LIB([lapack -lalapack -lcblas -latlas], [clapack_dgetrf],
	[AC_SUBST([LIBLAPACK_ATLAS], ["-llapack -lalapack -lf77blas"])
	AC_DEFINE([HAVE_LIBLAPACK_ATLAS], [1], [Use ATLAS LAPACK library])])])
AS_IF([test "x$with_atlas" != xno],
	[AC_CHECK_LIB([lapack -lf77blas -lgfortran -lcblas -latlas], [dsyev_],
	[AC_SUBST([LIBLAPACK_F77], ["-llapack -lf77blas -lgfortran"])
	AC_DEFINE([HAVE_LIBLAPACK_F77], [1], [Use LAPACK library])])])
AS_IF([test "x$with_atlas" != xno && test -z "$LIBATLAS"],
	[AC_MSG_FAILURE(
	[--with-atlas was given, but test for ATLAS failed])])
AS_IF([test "x$with_atlas" != xno && test -z "$LIBLAPACK_ATLAS"],
	[AC_MSG_FAILURE(
	[--with-atlas was given, but test for CLAPACK failed])])
AS_IF([test "x$with_atlas" != xno && test -z "$LIBLAPACK_F77"],
	[AC_MSG_FAILURE(
	[--with-atlas was given, but test for F77 LAPACK failed])])

AC_ARG_WITH([intel-mkl],
	[AC_HELP_STRING([--with-intel-mkl],
	[link against Intel Math Kernel Library])],
	[],
	[with_intel_mkl=no])

LIBMKL=
LIBMKL_LAPACK=
AS_IF([test "x$with_intel_mkl" != xno],
	[AC_CHECK_HEADER([mkl.h],
	AC_DEFINE([HAVE_MKL_H], [1], [Use mkl.h header]),
	[AC_MSG_FAILURE(
	[--with-intel-mkl was given, but test for mkl.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_intel_mkl" != xno],
	[AC_CHECK_HEADER([mkl_cblas.h],
	AC_DEFINE([HAVE_MKL_CBLAS_H], [1], [Use mkl_cblas.h header]),
	[AC_MSG_FAILURE(
	[--with-intel-mkl was given, but test for mkl_cblas.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_intel_mkl" != xno],
	[AC_CHECK_HEADER([mkl_lapack.h],
	AC_DEFINE([HAVE_MKL_LAPACK_H], [1], [Use mkl_lapack.h header]),
	[AC_MSG_FAILURE(
	[--with-intel-mkl was given, but test for mkl_lapack.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_intel_mkl" != xno],
	[AC_CHECK_LIB([mkl_intel_lp64 -lmkl_core -lmkl_sequential], [cblas_dgemm],
	[AC_SUBST([LIBMKL], ["-lmkl_intel_lp64 -lmkl_core -lmkl_sequential"])
	AC_DEFINE([HAVE_LIBMKL], [1], [Use Intel Math Kernel Library])])])
AS_IF([test "x$with_intel_mkl" != xno],
	[AC_CHECK_LIB([mkl_intel_lp64 -lmkl_core -lmkl_sequential], [dgetrf],
	[AC_SUBST([LIBMKL_LAPACK], ["-lmkl_intel_lp64 -lmkl_core -lmkl_sequential"])
	AC_DEFINE([HAVE_LIBMKL], [1], [Use Intel Math Kernel Library])])])
AS_IF([test "x$with_intel_mkl" != xno],
	[AC_CHECK_LIB([mkl_intel -lmkl_core -lmkl_sequential], [cblas_dgemm],
	[AC_SUBST([LIBMKL], ["-lmkl_intel -lmkl_core -lmkl_sequential"])
	AC_DEFINE([HAVE_LIBMKL], [1], [Use Intel Math Kernel Library])])])
AS_IF([test "x$with_intel_mkl" != xno],
	[AC_CHECK_LIB([mkl_intel -lmkl_core -lmkl_sequential], [dgetrf],
	[AC_SUBST([LIBMKL_LAPACK], ["-lmkl_intel -lmkl_core -lmkl_sequential"])
	AC_DEFINE([HAVE_LIBMKL], [1], [Use Intel Math Kernel Library])])])
AS_IF([test "x$with_intel_mkl" != xno &&
	( test -z "$LIBMKL" && test -z "$LIBMKL_LAPACK" )],
	[AC_MSG_FAILURE(
	[--with-intel-mkl was given, but test for Intel Math Kernel Library failed])])

AC_ARG_WITH([sunperf],
	[AC_HELP_STRING([--with-sunperf],
	[link against Sun Performance Library (only with Oracle Studio C compilers)])],
	[],
	[with_sunperf=no])

LIBSUNPERF=
AS_IF([test "x$with_sunperf" != xno],
	[AC_CHECK_HEADER([sunperf.h],
	AC_DEFINE([HAVE_SUNPERF_H], [1], [Use sunperf.h header]),
	[AC_MSG_FAILURE(
	[--with-sunperf was given, but test for sunperf.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_sunperf" != xno],
	[AC_CHECK_HEADER([cblas.h],
	AC_DEFINE([HAVE_CBLAS_H], [1], [Use cblas.h header]),
	[AC_MSG_FAILURE(
	[--with-sunperf was given, but test for cblas.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_sunperf" != xno],
	[CFLAGS="$CFLAGS -xlic_lib=sunperf"], [])
AS_IF([test "x$with_sunperf" != xno],
	[AC_CHECK_LIB([cblas], [cblas_dgemm],
	[AC_SUBST([LIBSUNPERF], ["-lcblas"])
	AC_DEFINE([HAVE_LIBSUNPERF], [1], [Use Sun Performance Library])],
	[AC_MSG_FAILURE(
	[--with-sunperf was given, but test for libsunperf failed])])])
AS_IF([test "x$with_sunperf" != xno],
	[AC_CHECK_LIB([cblas], [dgetrf],
	[AC_SUBST([LIBSUNPERF], ["-lcblas"])
	AC_DEFINE([HAVE_LIBSUNPERF], [1], [Use Sun Performance Library])],
	[AC_MSG_FAILURE(
	[--with-sunperf was given, but test for libsunperf failed])])])

AC_ARG_WITH([accelerate],
	[AC_HELP_STRING([--with-accelerate],
	[link against Apple Accelerate Framework (Mac OS X only)])],
	[],
	[with_accelerate=no])

LIBACCELERATE=
AS_IF([test "x$with_accelerate" != xno],
	[AC_CHECK_HEADER([clapack.h],
	AC_DEFINE([HAVE_CLAPACK_H], [1], [Use clapack.h header]),
	[AC_MSG_FAILURE(
	[--with-accelerate was given, but test for clapack.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_accelerate" != xno],
	[AC_CHECK_HEADER([cblas.h],
	AC_DEFINE([HAVE_CBLAS_H], [1], [Use cblas.h header]),
	[AC_MSG_FAILURE(
	[--with-accelerate was given, but test for cblas.h failed])],
	[AC_INCLUDES_DEFAULT])])
AS_IF([test "x$with_accelerate" != xno],
	[AC_CHECK_LIB([BLAS -lLAPACK], [cblas_dgemm],
	[AC_SUBST([LIBACCELERATE], ["-framework Accelerate"])
	AC_DEFINE([HAVE_LIBACCELERATE], [1], [Use Apple Accelerate Framework])],
	[AC_MSG_FAILURE(
	[--with-accelerate was given, but test for Accelerate Framework failed])])])
AS_IF([test "x$with_accelerate" != xno],
	[AC_CHECK_LIB([BLAS -lLAPACK], [dgetrf_],
	[AC_SUBST([LIBACCELERATE], ["-framework Accelerate"])
	AC_DEFINE([HAVE_LIBACCELERATE], [1], [Use Apple Accelerate Framework])],
	[AC_MSG_FAILURE(
	[--with-accelerate was given, but test for Accelerate Framework failed])])])

AS_IF([test "x$with_atlas" = xno] \
	&& [test "x$with_intel_mkl" = xno] \
	&& [test "x$with_sunperf" = xno] \
	&& [test "x$with_accelerate" = xno],
	[AC_MSG_FAILURE( \
	[please indicate one of the following: --with-atlas, --with-intel-mkl, \
--with-sunperf, --with-accelerate])])

# Checks for buggy rint()
# see https://bugzilla.redhat.com/show_bug.cgi?id=479025
# http://bugs.centos.org/view.php?id=3326
# http://sourceware.org/bugzilla/show_bug.cgi?id=5857

cat << eof > src/test.c
#include <math.h>

int main(int argc, char **argv)
{
	return (((int)rint(289440.75) == 289440) ? 0 : 1);
}
eof
AC_MSG_CHECKING(if rint() is buggy)
BUGGY_RINT=no
$CC $CFLAGS src/test.c -o src/test -lm 2> /dev/null
AS_IF([src/test],
	[AC_SUBST([BUGGY_RINT], ["yes"])], [AC_MSG_RESULT([no])])
rm -f src/test.c src/test
AS_IF([test "x$BUGGY_RINT" = xyes],
	[AC_MSG_RESULT(yes)])
AS_IF([test "x$BUGGY_RINT" = xyes],
	[AC_MSG_NOTICE([A safe_rint() replacement will be used instead of buggy rint()])])
AS_IF([test "x$BUGGY_RINT" = xyes],
	[AC_DEFINE([HAVE_SAFE_RINT], [1], [Use safe_rint() replacement for buggy rint()])])

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile])
AC_OUTPUT
