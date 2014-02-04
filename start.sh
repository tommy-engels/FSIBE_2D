export LC_ALL=C
unset MAKEFLAGS
module load intel/11.1.080
OAR_FILE="meso_OMP.oar"
#OAR_FILE="job.oar"


#------------------------------------------------------
# COLOR, baby!
#------------------------------------------------------
# Reset
Color_Off='\e[0m'       # Text Reset

# Regular Colors
Black='\e[0;30m'        # Black
Red='\e[0;31m'          # Red
Green='\e[0;32m'        # Green
Yellow='\e[0;33m'       # Yellow
Blue='\e[0;34m'         # Blue
Purple='\e[0;35m'       # Purple
Cyan='\e[0;36m'         # Cyan
White='\e[0;37m'        # White

# Bold
BBlack='\e[1;30m'       # Black
BRed='\e[1;31m'         # Red
BGreen='\e[1;32m'       # Green
BYellow='\e[1;33m'      # Yellow
BBlue='\e[1;34m'        # Blue
BPurple='\e[1;35m'      # Purple
BCyan='\e[1;36m'        # Cyan
BWhite='\e[1;37m'       # White

#------------------------------------------------------
# check call ok
#------------------------------------------------------
if [ "$1" == "" ]
then
  echo -e ${BCyan}"make or not?"${Color_Off}
  exit 1
fi

if [ "$2" == "" ]
then
  echo -e ${BCyan}"directory?"${Color_Off}
  exit 1
fi

#------------------------------------------------------
# check if directory exists.
#------------------------------------------------------
if [ $2 != '' ]
then
  if [ ! -d "$2" ]
    then
    echo -e ${BRed}"did not find directory"${Color_Off}
    exit 1
  fi  
fi


#------------------------------------------------------
# make if you like
#------------------------------------------------------
if [ $1 == "make" ]
then 
  make meso
fi

# success?
if [ ! -f "dns.out" ] 
then 
  echo -e ${BRed}"MAKE FAILED"${Color_Off}
  exit 1
fi

akt=$PWD #get parent directory





#------------------------------------------------------
# clean directory
#------------------------------------------------------
if [ "$3" == "yes" ] # clean directory if desired
then 
  cd "$2"
  rm -Rf vor/
  rm -Rf vor2/
  rm -Rf press/
  rm -Rf press2/
  rm -Rf fields/
  rm -Rf mask/
  rm -Rf mvt/
  rm -Rf fix/
  mkdir tmp
  mv *.oar tmp/
  mv PARAMS.m tmp/
  mv epsilons tmp/
  mv *.in tmp/
  mv *.inicond tmp/
  mv *.sh tmp/
  rm *
  mv tmp/* ./
  rm -R tmp/
  cd $akt #go back to parent directory
fi





echo -e ${BCyan}"------------------------------------"${Color_Off}
echo -e ${BGreen}$2${Color_Off}
# copy exceutable file 
cp dns.out "$2"
archiv="$(date +'%Y.%m.%d_%H.%M.%S').code.tar.gz"
echo ${archiv}
tar czf "${archiv}" *.f90 makefile PARAMS.m
mv ${archiv} $2
cd "$2" 
#--now we're in the target directory
echo -e ${BCyan}"------------------------------------"${Color_Off}




if [ "$4" == "sure" ]
then
    #--------------------------------------
    # Start directly without asking
    #--------------------------------------
    echo -e ${Yellow}
    cat ${OAR_FILE}
    echo -e ${Cyan}
    cat PARAMS.m
    echo -e ${Color_Off}
    
    oarsub -S  ./${OAR_FILE}
    echo -e ${BGreen}"*** RUNNING ***"${Color_Off}
else
    #--------------------------------------
    # Start, but interactively
    #--------------------------------------
    start='k'
    while [ $start != "n" ]
    do 
      echo -e ${Yellow}
      cat ${OAR_FILE}
      echo -e ${Cyan}
      cat PARAMS.m
      echo -e ${Color_Off}
      
      echo -n -e ${BGreen}"Start ? (y for start, P to edit PARAMS.m, D to edit dns.pbs, n to abort "${Color_Off}
      read start
      
      if [ $start == "y" ]
      then 
	start='n'        
	oarsub -S  ./${OAR_FILE}
	echo -e ${BGreen}"*** RUNNING ***"${Color_Off}
      fi
      if [ $start == "P" ]
      then
	nano PARAMS.m
      fi
      if [ $start == "D" ]
      then
	nano ${OAR_FILE}
      fi      
    done     
fi





  
