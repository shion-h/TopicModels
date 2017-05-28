#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(TopicModels, 1.0, [dummy@xxx.yyy.zz])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/AuthorFileParser.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lboost_program_options':
AC_CHECK_LIB([boost_program_options], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
