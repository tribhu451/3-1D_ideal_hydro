#!/bin/bash



function hadronisation()
{
rm -rf therminator/fomodel/lhyquid3d/*
echo "Now no file is there inside therminator/fomodel/lhyquid3d/ "
cp -r hydro_output/*.xml therminator/fomodel/lhyquid3d/hypersurface.xml
echo "moving the .xml file from hydro_output to therminator/fomodel/lhyquid3d/ "
cd therminator
echo "entered into therminator"
rm -rf ./events/* 
echo "removed all the contents of folder therminator/events/ (the folder where output will be stored)."
nohup ./therm2_events > ../logFiles/therminator.log &
echo "job submitted in therminator. "
cd ..
echo "exiting from therminator ..."
}






# The flag file
file0='hydro_complete_flag.txt'

# executable file of hydro
file1='3Dhydro'

# [step 1] delete the flag file (if exist ...)
if [ -f $file0 ]; then
     rm -rf $file0
     echo "initially there was a flag file, now it's deleted ... "
else
     echo "there was no flag file initially"
fi



# [step 2] compilaton
make clean
make cfiles
make
sleep 3


# [step 3] submitting the job
if [ -f $file1 ]; then
    nohup ./$file1 > ./logFiles/hydro_run.log &
    echo "job submitted ..."
    # [step 3.1] when the job will be completed
    count=0
    while [ $count -lt 36000 ]
    do
       if [ -f $file0 ]; then
           echo "hydro job completed in $count sec."
           echo "Now going to perform hadronisation (THERMINATOR)"
           hadronisation 
           break
       else
           ((count=count+210))
           echo "wait more 3 min. 30 sec...."
           sleep 210
       fi
     done

else
  # [step 3.2] if compilation error occurs...
  echo "compilation error ..."
  exit

fi


