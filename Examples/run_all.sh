#!/usr/bin/env bash
SLR=`realpath -e ../bin/Slr`
NPARALLEL=2

ls -d * | xargs -n 1 -I %% -P ${NPARALLEL} sh -c "
        echo '* Running %%'
	cd %%
	if [ ! -e 'SKIP_RUN' ]
	then 
	    OPENBLAS_NUM_THREADS=1 stdbuf -o 0 ${SLR} slr.ctl > run.log
	fi
"
