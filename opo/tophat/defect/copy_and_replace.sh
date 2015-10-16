#!/bin/bash

sed -e 's/#.*$//' -e '/^\s*$/d' ${1} | awk '{print $'${3}'}' | awk '{if(NR % '${2}' == 0){print;}else{printf "%s ",$0}}'| tac  > ${1}.copy 
