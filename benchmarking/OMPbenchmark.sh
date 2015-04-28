#!/bin/bash
#PBS -l procs=16,tpn=16,mem=46gb,walltime=2:00:00
#PBS -q cpsc424

# For PBS scripts
if [ -n "$PBS_O_WORKDIR" ]; then cd $PBS_O_WORKDIR; fi

BASEDIR=~/scratch/SupportEnumeration
GAMESDIR=$BASEDIR/games
EXEC=$BASEDIR/OMPsupEnum
LOGFILE=OMPdata.log
TESTFILES=testFiles.txt

trap "echo 'stopping OMP script'; exit" SIGINT SIGTERM

# Takes 1 arg, and takes contents of arg and writes to log file
log () {
    echo $1 >> $LOGFILE
}

echo $(date)": start benchmarking for OMPsupEnum using games in file $TESTFILES"

log $(date)": start benchmarking for OMPsupEnum using games in file $TESTFILES"
for test_file in `cat $TESTFILES`; do      # Exec trials for each test file
    log "$test_file: "
    for threads in 2 4 8 16; do
        log "For $threads threads:"
        for i in {1..2}; do
            OUT=$($EXEC $GAMESDIR/$test_file $threads)
            log "Trial $i: $OUT"                        # Record return 
                                                        # val (timing) 
        done
    done
done

