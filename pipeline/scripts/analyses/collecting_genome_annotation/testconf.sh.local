#!/usr/bin/sh
# Test if softwraes are available

# testing test

which programe_bidon
status=$?
echo "Test bidon : $status"
if  [ $status -eq 1 ]
then
echo "Ok"
else
echo "Error."
exit
fi



echo "MISSING SOFTWARE " > missing_soft 
echo "================ " >> missing_soft 



# Executables
# ===========

for soft in "esearch" "gffread" "wget" "gzip" 
do

which $soft
status=$?
echo "Test $soft : $status"
if  [ $status -eq 0 ]
then
echo "Ok"
else
echo "Error, $soft is missing"
echo $soft >> missing_soft 
fi
done







 


