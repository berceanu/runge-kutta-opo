#!/bin/bash

sed -e 's/#.*$//' -e '/^\s*$/d' ${1} | awk '{print $2}' | tac > ${1}.copy
