#!/bin/sh
unset LANG
ulimit -s unlimited

IFORT="ifort"
COMPILO=$IFORT


NETCDFLIB=/opt/netcdf/$COMPILO/lib
NETCDFINC=/opt/netcdf/$COMPILO/include


NAME=$1
instant=$2
nol=$3

FC4_BUG=-no-ipo
prog="A-GCM"


F90FLAGS1="-I${HDFINC} -L${HDFLIB} -ljpeg -lz -I${NETCDFINC} -L${NETCDFLIB} -lnetcdf" 
CC="cc -DLITTLE"
CFLAGS="-c -g -DUNDERSCORE -DLITTLE"


VERSION=`echo $NAME | cut -c7-10`
 
if [ -z $instant ]; then

sed s/Prog_version/$NAME/g $NAME.f90 > $NAME.tmp.f90 

else

if (( $nol == 1)); then
nol2="nol"
sed s/Prog_version/$NAME/g $NAME.f90 | sed s/Num_version/$VERSION$nol2/ | sed s/instant_switch/$instant/g | sed s/nol_switch/$nol/ > $NAME.tmp.f90 

else
nol2=""
sed s/Prog_version/$NAME/g $NAME.f90 | sed s/Num_version/$VERSION$nol2/ | sed s/instant_switch/$instant/g | sed s/nol_switch/$nol/ > $NAME.tmp.f90 

fi

fi


${COMPILO} $NAME.tmp.f90 $F90FLAGS1 -traceback -o ${NAME}.e


rm -f $NAME.tmp.f90
