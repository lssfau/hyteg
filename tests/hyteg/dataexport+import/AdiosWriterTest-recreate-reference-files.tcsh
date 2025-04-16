#!/bin/tcsh

echo "Going to generate BP files:"

mpirun -n 1 ./AdiosWriterTest
mpirun -n 3 ./AdiosWriterTest

cd AdiosWriterTest-Output

foreach BP (`/bin/ls -d *.bp` )
    set output=`basename $BP .bp`.txt
    echo "---> Creating "$output
    bpls -la $BP > $output
end

cd ..

echo " "
echo " ... done "
echo " "
echo "Please check txt-files and copy them into"
echo "<src-dir>/hyteg/dataexport+import/AdiosWriterTest-References"
echo " "
echo "Take care to check for / remove blank lines at the file ends!"
echo " "
