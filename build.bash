#!/bin/bash -e

mydate=`date "+%D %H:%M:%S"`

echo "static const char *build_date = \"$mydate\";" \
  > build_date.txt

txt2cpp.py help.txt \
  > usage.h

vcxproj2makefile.py

make
