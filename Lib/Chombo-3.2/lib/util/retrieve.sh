#!/bin/sh

if [ -x "`which curl`" ] 
  then curl $1 -o $2
elif [ -x "`which wget`" ]
  then wget -nc $1 -O $2
else
  echo "Could not find curl or wget!"
  exit 1
fi
