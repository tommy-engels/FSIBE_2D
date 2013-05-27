echo "====================================================="
echo "==                                                 =="
echo "==                                                 =="
echo "==                                                 =="
echo "==                                                 =="
echo "==                                                 =="
echo "==                                                 =="
echo "====================================================="


export LC_ALL=C

if [ $2 != '' ]
then
    if [ -d "$2" ]
    then
      akt=$PWD #get parent directory
      if [ $1 == "make" ]
      then 
        make duke
      fi
      if [ -f dns.out ] 
      then         
        if [ $3 == "yes" ] # clean directory if desired
        then 
          cd "$2"
          rm -R vor/
          rm -R vor2/
          rm -R press/
          rm -R press2/
          rm -R fields/
          rm -R mvt/
          rm -R fix/
          mkdir tmp
          mv duke.oar tmp/
          mv PARAMS.m tmp/
          mv epsilons tmp/
	  mv *.in tmp/
	  mv *.sh tmp/
          rm *
          mv tmp/* ./
          rm -R tmp/
          cd $akt #go back to parent directory
        fi  # end clean directory
          echo -e "\033[49;5;32m ========================path: $2 \033[0m"
          cp dns.out "$2"
          archiv="$(date +'%Y.%m.%d_%H.%M.%S').code.tar.gz"
          echo $archiv
          tar czf "$archiv" *.f90 makefile PARAMS.m
          mv "$archiv" "$2"
          cd "$2"
          echo -e "\033[49;5;32m ================================================================ \033[0m"
          if [ "$4" == "sure" ]
          then
              cat duke.oar
              cat PARAMS.m
              oarsub -S  ./duke.oar
              echo -e "\033[49;5;32m *** RUNNING *** \033[0m"
          else
              more dns.pbs PARAMS.m
              start='k'
              while [ $start != "n" ]
              do 
                echo -n -e "\033[49;5;32m Start ? (y for start, P to edit PARAMS.m, D to edit dns.pbs, n to abort \033[0m"
                read start
                if [ $start == "y" ]
                then 
                  start='n'        
                  oarsub -S  ./duke.oar
                fi
                if [ $start == "P" ]
                then
                  nano PARAMS.m
                fi
                if [ $start == "D" ]
                then
                  nano duke.oar
                fi      
              done    
              cd ..
          fi

      else
	    echo -e "\033[49;5;32m make a échoué \033[0m"
      fi
    else
    echo -e "\033[49;5;32m directory inexistent \033[0m"
    fi
else
    echo -e "\033[49;5;32m no directory specified \033[0m"
fi
