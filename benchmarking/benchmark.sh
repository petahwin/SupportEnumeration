#!/bin/bash

trap "echo 'stopping general script'; exit" SIGINT SIGTERM

./Serialbenchmark.sh
./OMPbenchmark.sh

