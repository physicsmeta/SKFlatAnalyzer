#!/bin/bash

cd /data6/Users/jihkim/SKFlatRunlog/2019_10_16_205309__168954__ChargeFlipValidation__Year2016__CFrate__TAMSA1/DYJets 

myvar=`grep -rF [SKFlatNtuple::Loop\ RUNNING]\ 0`
myvar1=`grep -r electrons\ :\ 1`;
myvar2=`grep -r electrons\ :\ 2`;
myvar3=`grep -r electrons\ :\ 3`;

count=0;
count1=0;
count2=0;
count3=0;

for var in $myvar; do
  if [[ $var = 0/* ]]
    then count_tmp=${var#0/}
        #echo count_tmp : $count_tmp
         count=`expr $count + $count_tmp`
        #echo count : $count
  fi
done

for var1 in $myvar1; do
  if [[ $var1 = job* ]]
    then ((count1++))
  fi
done

for var2 in $myvar2; do
  if [[ $var2 = job* ]]
    then ((count2++))
  fi
done

for var3 in $myvar3; do
  if [[ $var3 = job* ]]
    then ((count3++))
  fi
done

RED='\033[0;31m'
LRED='\033[0;91m'
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
NC='\033[0m'

echo -e ${RED}\[ChargeFlipValidation] ${LRED}CFrate w/ no requirement on the \# of electrons${NC}
echo -e ${YELLOW}total events${NC} : $count
echo -e ${BLUE}events with reco electrons 1${NC} : $count1 
echo -e ${BLUE}events with reco electrons 2${NC} : $count2 
echo -e ${BLUE}events with reco electrons 3${NC} : $count3 



