#!/bin/bash
#num=$(awk 'BEGIN{for(i=1;i<=3;i+=0.2)print i}')
#num=(1.425 2.525 7.25,7.5,7.75)
#for j in $num;
#num=$(awk 'BEGIN{for(i=1;i<=10;i+=0.1)print i}')
#for j in 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.4 3.6 3.8 4.0 4.2 4.4 4.6 4.8 5.0 5.2 5.4 5.6 5.8 6.0 6.5 7.0 7.5 8.0 9.0 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
for j in 5.0 5.2 5.4 5.6 5.8 6.0 6.5 7.0 7.5 8.0 9.0 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 30 40 50 60 80 100 120                            
#for j in  60 80

#for j in 50 60 80
#for j in 0.4 0.7 1.75
#for j in 3.2 3.4 3.6 3.8 4.2 4.4 4.6 4.8
#for j in 5.2 5.4 5.6 5.8 #4.2 4.4 4.6 4.8
#for j in $(seq 11 20);
#for j in 1.85 1.9 1.95 2.05 2.1 2.15
#for j in  2.85 2.9 2.95 23
do
    echo "Creating directory..."
    m_newrun template $j
    cd $j
    sed -i -e "s/initial_mass =/initial_mass = $j /" inlist_grid
    echo "Running Model..."
    ./rn > rn.out
    cd ..
done
