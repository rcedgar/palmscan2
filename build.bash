#!/bin/bash -e

mydate=`date "+%D %H:%M:%S"`
githash=`git rev-parse --short HEAD`

echo "static const char *build_info = \"$mydate:$githash\";" \
  > build_info.txt

txt2cpp.py usage.txt \
  > usage.h

vcxproj2makefile.py

make
