#!/bin/bash
#PBS -l procs=1,tpn=1,mem=46gb,walltime=2:00:00
#PBS -q cpsc424

# For PBS scripts
if [ -n "$PBS_O_WORKDIR" ]; then cd $PBS_O_WORKDIR; fi

BASEDIR=~/scratch/SupportEnumeration
GAMESDIR=$BASEDIR/games
EXEC=$BASEDIR/SerialsupEnum
LOGFILE=Serialdata.log
TESTFILES=testFiles.txt

trap "echo 'stopping serial script'; exit" SIGINT SIGTERM

# Takes 1 arg, and takes contents of arg and writes to log file
log () {
    echo $1 >> $LOGFILE
}

echo "$(date): start benchmarking for SerialsupEnum using games in file $TESTFILES"

log "$(date): start benchmarking for SerialsupEnum using games in file $TESTFILES"
for test_file in `cat $TESTFILES`; do      # Exec trials for each test file
    log "$test_file: "
    for i in {1..2}; do
        OUT=$($EXEC $GAMESDIR/$test_file)
        log "Trial $i: $OUT"                # Record return val (timing) 
                                            # of trial of execution
    done
done

