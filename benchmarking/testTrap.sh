#!/bin/bash

trap "echo 'hello world'; exit" SIGINT SIGTERM

for i in {1..10}; do
    echo "sleeping $i..."
    sleep 5
done

