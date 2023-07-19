#!/bin/bash

for i in 2 3 4 5 6 7 8; do
    cp d1/run.sh d$i;
    perl -i -p -e "s/drone=1/drone=${i}/g" d$i/run.sh
    echo "updated d$i/run.sh"
done

	 
