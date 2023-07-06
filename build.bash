#!/bin/bash -e

mydate=`date "+%D %H:%M:%S"`

echo "static const char *build_date = \"$mydate\";" \
  > build_date.txt

vcxproj2makefile.py

make
