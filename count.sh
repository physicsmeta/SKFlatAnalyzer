#!/bin/bash

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
         echo count_tmp : $count_tmp
         count=`expr $count + $count_tmp`
         echo count : $count
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

echo total events : $count
echo events with reco electrons 1 : $count1
echo events with reco electrons 2 : $count2
echo events with reco electrons 3 : $count3



