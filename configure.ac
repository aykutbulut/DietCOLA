AC_PREREQ([2.69])
AC_INIT([dietcola], [0.0], [aykut@lehigh.edu])
AC_CONFIG_SRCDIR([config.h.in])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_HEADERS([config.h])


# for creating library
AC_PROG_LIBTOOL

# Checks for programs.
AC_PROG_CXX

# Installs to ./bin, changes default behavior
AC_PREFIX_DEFAULT([`pwd`])

# Checks for header files.
AC_CHECK_HEADERS([math.h])

# Checks for typedefs, str
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Check dependent libraries
#AC_CHECK_LIB(library, function, [action-if-found], [action-if-not-found], [other-libraries])
AC_CHECK_LIB([m], [sqrt])
PKG_CHECK_MODULES([COLA], [cola >= 0.0])
PKG_CHECK_MODULES([OSIMOSEK], [osimosek >= 0.0])
PKG_CHECK_MODULES([OSICONIC], [osiconic >= 0.0])
PKG_CHECK_MODULES([ALPS], [alps >= 1.4.8])
PKG_CHECK_MODULES([CGL], [cgl >= 0.58.5])

AC_CONFIG_FILES([Makefile
	src/Makefile
	dietcola.pc
	dietcola-uninstalled.pc])

AC_OUTPUT
