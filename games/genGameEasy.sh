#!/bin/bash

cmdlineargs="$*"

# Generate the game
java -jar gamut.jar $cmdlineargs

# Continue only if game was successfully generated
if [ "$?" -eq 0 ]; then

# Find the file name of the game, pass it into  conversion script
./convEasyParse.py `find -type f -printf '%TY-%Tm-%Td %TT %p\n' | sort -r | \
head -n 1 | awk '{print $3}'`

else

echo "Script failed at gamut.jar call"

fi
