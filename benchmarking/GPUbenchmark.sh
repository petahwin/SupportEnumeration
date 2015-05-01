#!/bin/bash
#PBS -l procs=1,tpn=1,mem=46gb,walltime=1:00:00
#PBS -q phitest

# For PBS scripts
if [ -n "$PBS_O_WORKDIR" ]; then cd $PBS_O_WORKDIR; fi

BASEDIR=~/scratch/SupportEnumeration
GAMESDIR=$BASEDIR/games
EXEC=$BASEDIR/GPUsupEnum
LOGFILE=GPUdata.log
TESTFILES=testFiles.txt

trap "echo 'stopping GPU script'; exit" SIGINT SIGTERM

# Takes 1 arg, and takes contents of arg and writes to log file
log () {
    echo $1 >> $LOGFILE
}

echo "$(date): start benchmarking for GPUsupEnum using games in file $TESTFILES"

log "$(date): start benchmarking for GPUsupEnum using games in file $TESTFILES"
for test_file in `cat $TESTFILES`; do      # Exec trials for each test file
    log "$test_file: "
    for i in {1..2}; do
        OUT=`{ time $EXEC $GAMESDIR/$test_file ; } 2>&1`
        log "Trial $i: $OUT"                        # Record return 
                                                    # val (timing) 
    done
done

