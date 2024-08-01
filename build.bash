#!/bin/bash -e

mydate=`date "+%D %H:%M:%S"`
githash=`git rev-parse --short HEAD`

echo "static const char *build_info = \"$mydate:$githash\";" \
  > build_info.txt

python3 $py/txt2cpp.py usage.txt \
  > usage.h

python3 $py/vcxproj2makefile.py

make
