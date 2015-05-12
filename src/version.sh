#!/bin/sh

ESCAPESCRIPT='s/\([\\"]\)/\\\\\1/g'

####################
### Gather flex info
if test -z "${LEX}"; then
    LEX="flex"
fi
LEX_PATH=`which "${LEX}" 2>/dev/null`
if test -z "${LEX_PATH}"; then
  LEX="(missing)"
  LEX_PATH="(flex not found)"
elif test -x "${LEX_PATH}"; then
    LEX_VERSION=`${LEX} --version`
else
    LEX_LEXER_PREBUILT=`dirname $0`/mdllex.flex.c
    if test -r "${LEX_LEXER_PREBUILT}"; then
        LEX_VERSION="not installed - using pregenerated lexer"
    else
        LEX_VERSION="not installed"
    fi
fi
LFLAGS="-Crema ${LFLAGS}"
LFLAGS=`echo $LFLAGS | sed "${ESCAPESCRIPT}"`

####################
### Gather bison info
if test -z "${YACC}"; then
  YACC="bison"
elif test "${YACC}" = "bison -y"; then
  YACC="bison"
fi
YACC_PATH=`which "${YACC}" 2>/dev/null`
if test -z "${YACC_PATH}"; then
  YACC="(missing)"
  YACC_PATH="(bison not found)"
elif test -x "${YACC_PATH}"; then
    YACC_VERSION=`${YACC} --version | head -1`
else
    YACC_PARSER_PREBUILT=`dirname $0`/mdlparse.bison.c
    if test -r "${YACC_PARSER_PREBUILT}"; then
        YACC_VERSION=`head -1 ${YACC_PARSER_PREBUILT} | sed 's/^.*GNU Bison \(.*\)\. .*/not installed - using parser pregenerated by GNU Bison \1/'`
    else
        YACC_VERSION="not installed - build should fail"
    fi
fi
if test -z "${YFLAGS}"; then
    YFLAGS=""
fi
YFLAGS=`echo $YFLAGS | sed "${ESCAPESCRIPT}"`

####################
### Gather cc info
if test -z "${CC}"; then
    CC="cc"
fi
CC_PATH=`which ${CC} 2>/dev/null`
if test -x "${CC_PATH}"; then
    CC_VERSION=`${CC} --version | head -1`
fi
if test -z "${CFLAGS}"; then
    CFLAGS=""
fi
CFLAGS=`echo $CFLAGS | sed "${ESCAPESCRIPT}"`
if test -z "${LD}"; then
    LD="ld"
fi

####################
### Gather ld info
LD_PATH=`which ${LD} 2>/dev/null`
if test -x "${LD_PATH}"; then
    LD_VERSION=`${LD} --version | head -1`
fi
if test -z "${LDFLAGS}"; then
    LDFLAGS=""
fi
LDFLAGS=`echo $LDFLAGS | sed "${ESCAPESCRIPT}"`

####################
### Gather buildhost info
BUILD_USER=`id -un`
BUILD_HOST=`hostname -f`
SRC_DIR=`dirname $0`
SRC_DIR=`sh -c "cd ${SRC_DIR} && pwd"`
BUILD_DIR=`pwd`
BUILD_DATE=`date`
BUILD_UNAME=`uname -a`

####################
### Gather program version info
version_path=`dirname $0`/version.txt
MCELL_VERSION=`cat "${version_path}"`
MCELL_HAVE_REVISION_INFO=0
GIT=`which git 2>/dev/null`
if test -x "${GIT}"; then
  if "${GIT}" describe --all --abbrev=4 HEAD >/dev/null 2>&1; then
    MCELL_HAVE_REVISION_INFO=1
  fi
fi

echo "/****************************************************************"
echo " * This file is automatically generated.  Do not edit it."
echo " * version.h updated at ${BUILD_DATE} on ${BUILD_HOST} by ${BUILD_USER}"
echo " ****************************************************************/"
echo
echo "/* Program version info */"
echo "#define MCELL_VERSION \"${MCELL_VERSION}\""
if test "$MCELL_HAVE_REVISION_INFO" = "1"; then
    echo "#define MCELL_REVISION $(git show -s --format=\"%h\")"
    echo "#define MCELL_REVISION_DATE $(git show -s --format=\"%aD\")"
    echo "#define MCELL_REVISION_COMMITTED 1"
    echo "#define MCELL_REVISION_BRANCH \"$(git rev-parse --abbrev-ref HEAD)\""
else
    echo "#define MCELL_REVISION \"\""
    echo "#define MCELL_REVISION_DATE \"\""
    echo "#define MCELL_REVISION_COMMITTED 0"
    echo "#define MCELL_REVISION_BRANCH \"not in VCS\""
fi
echo
echo "/* Build info */"
echo "#define MCELL_BUILDDATE  \"${BUILD_DATE}\""
echo "#define MCELL_BUILDUSER  \"${BUILD_USER}\""
echo "#define MCELL_BUILDHOST  \"${BUILD_HOST}\""
echo "#define MCELL_SRCDIR     \"${SRC_DIR}\""
echo "#define MCELL_BUILDDIR   \"${BUILD_DIR}\""
echo "#define MCELL_BUILDUNAME \"${BUILD_UNAME}\""
echo
echo "/* Tool identity and version info */"
echo "#define MCELL_FLEX \"${LEX}\""
echo "#define MCELL_FLEX_PATH \"${LEX_PATH}\""
echo "#define MCELL_FLEX_VERSION \"${LEX_VERSION}\""
echo "#define MCELL_FLEX_FLAGS \"${LFLAGS}\""
echo "#define MCELL_BISON \"${YACC}\""
echo "#define MCELL_BISON_PATH \"${YACC_PATH}\""
echo "#define MCELL_BISON_VERSION \"${YACC_VERSION}\""
echo "#define MCELL_BISON_FLAGS \"${YFLAGS}\""
echo "#define MCELL_CC \"${CC}\""
echo "#define MCELL_CC_PATH \"${CC_PATH}\""
echo "#define MCELL_CC_VERSION \"${CC_VERSION}\""
echo "#define MCELL_LD \"${LD}\""
echo "#define MCELL_LD_PATH \"${LD_PATH}\""
echo "#define MCELL_LD_VERSION \"${LD_VERSION}\""
echo
echo "/* Build options */"
echo "#define MCELL_LFLAGS \"${LFLAGS}\""
echo "#define MCELL_YFLAGS \"${YFLAGS}\""
echo "#define MCELL_CFLAGS \"${CFLAGS}\""
echo "#define MCELL_LDFLAGS \"${LDFLAGS}\""
