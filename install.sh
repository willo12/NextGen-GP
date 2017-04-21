#!/usr/bin/env bash

if ! [ -d ${HOME}/bin ]; then
  mkdir ${HOME}/bin
fi

dir=$(pwd)
echo 'Creating ~/bin sym links to files in ' ${dir}

for e in disp driver find_elite runtree;do 
  echo "${HOME}/bin/${e} --> ${dir}/${e}.py";
  ln -s ${dir}/${e}.py ${HOME}/bin/${e};
done
