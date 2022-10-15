#!/bin/bash

MACHINE_ID=`uname -n`

if [ ${MACHINE_ID:0:3} == "MBP" ]
then

   git -C /Users/21b/Desktop/KORC rev-parse HEAD > git_hash.txt

   git -C /Users/21b/Desktop/KORC diff HEAD > git_diff.txt
   
elif [ ${MACHINE_ID:0:3} == "nid" ] || [ ${MACHINE_ID:0:4} == "cori" ]
then

   git -C /global/cfs/cdirs/m3236/build_unstable/KORC rev-parse HEAD > git_hash.txt

   git -C /global/cfs/cdirs/m3236/build_unstable/KORC diff HEAD > git_diff.txt
    
else
# MACHINE_ID is new and unknown. Inform the user how to add support for this new machine.
    echo $MACHINE_ID not suported by this script.
    echo To support this machine, add a new elif statement of the form
    echo
    echo elif [ \$MACHINE_ID == \"$MACHINE_ID\" ]
fi

exit $?
