#!/bin/bash

list="
myfunc
mypro
pssm2Qij
pssm2modm
rootname
search_new
subset
"

i=0
for rtname in $list ; do 
    ((i++))
    echo "=== $i ==="
    diffcontent=` diff $rtname.cpp ../../frag1d/src/$rtname.cpp` 
    if [ "$diffcontent" ]; then
        echo "$rtname.cpp differ"
        echo "run vimdiff  $rtname.cpp ../../frag1d/src/$rtname.cpp"
    fi
done
