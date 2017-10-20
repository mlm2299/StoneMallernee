#!/bin/bash

rm -r batch_infiles
mkdir batch_infiles

python main_big_in.py
cp ./big.in batch_infiles/big.in
rm big.in

for ((i=0; i<32; i++));
do
###{
###while ([ ! -f /batch_infiles/small.in${i} ]);
###do
###echo "making small.in${i}"
python main_small_in.py
cp ./small.in batch_infiles/small.in${i} 
rm small.in
###done
###}
done

wait
