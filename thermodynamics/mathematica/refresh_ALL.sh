#!/bin/sh
[ ! -e ALL.math ] && touch ALL.math
if [ -w ALL.math ]; then
  rm -f ALL.math; l=`ls -1 *.math`; dir=`pwd`
  echo allModulesLoaded=True > ALL.math;
  for i in $l; do echo "<<$dir/$i"; done >> ALL.math
else
  echo; echo -e "\033[31m\033[1mWARNING: no write permission for refresh_ALL.math\033[0m" 2>&1
fi
