#!/bin/bash

#  $Id: envGMC.sh 268 2013-10-31 17:21:49Z tassielli $
#  $Author:  $
#  $Revision:  $

# base PRJ
if [ -v PRJBASE ];
then
	echo "PRJ BASE already set to "${PRJBASE}
else
	export PRJBASE="/mnt/sndhd/sw/FIRB-sw"
fi

#gcc xml

#geant4
cd $G4INSTBASE/share/Geant4-10.1.3/geant4make
source geant4make.sh
cd -

#set the workdir for g4GMC
#export homedir=`pwd`
export G4WORKDIR=${PRJBASE}/simulation/g4GMC
export LD_LIBRARY_PATH=${G4WORKDIR}/lib:${PRJBASE}/simulation/ConfigReader/lib:${LD_LIBRARY_PATH}

#rome
if [ -v ROMESYS ]
then
	echo "ROME already set to "${ROMESYS}
else
 export ROMESYS=/mnt/sndhd/sw/meg/rome
 export PATH=$ROMESYS/bin:${PATH}
 export LIBROME=yes
fi

#root
#export ROOTSYS=/pro/root_v5_34_30
#export PATH=$ROOTSYS/bin:${PATH}
#export LD_LIBRARY_PATH=${ROOTSYS}/lib:${LD_LIBRARY_PATH}

#midas
#export MIDASSYS=/mnt/sndhd/sw/meg/midas
#export MIDAS_DIR=/mnt/sndhd/sw/meg/midas
#export MIDAS_EXPTAB=$MIDAS_DIR/exptab
#export MIDAS_EXPT_NAME=FIRB_DAQ
#export PATH=${MIDASSYS}/linux/bin:$PATH
#export LD_LIBRARY_PATH=${MIDASSYS}/linux/lib:${LD_LIBRARY_PATH}

#GMC
if [ -v GMCDIR ]
then
	echo "GMCDIR already set to "${GMCDIR}
else
	export GMCDIR=${PRJBASE}/analyzer/GMC
	export LD_LIBRARY_PATH=${GMCDIR}/obj:${LD_LIBRARY_PATH}
fi
