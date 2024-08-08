#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
# Ask via console the start, end and step time. 
# The .0000 defines the input age format of the files
echo 'Enter the initial time with four decimal digits in Ma (e.g. 0.0000)'
read timeStart
echo 'Enter the end time with four decimal digits in Ma (e.g. 1000.0000)'
read timeEnd
echo 'Enter the step in Ma (e.g. 1)'
read step
# create a temporary file that has the header needed for the final result
cat << EOF > tmpHeader.txt
{
    "type": "AllAgesNdPolygons",
    "ages": [
EOF
# loop over all the files and concatenate them. The files are located in the "topologyFiles" folder and with a name of "topology_AgeMa.geojson", example: topologyFiles/topology_1.0000Ma.geojson.
for i in $(seq $timeStart $step $timeEnd) 
do  
	cat "topologyFiles/topology_"$i"Ma.geojson"
 # when you reach to the end tipe, close the parenthesis and comas
	if [[ "$i" != "$timeEnd" ]]; then 
		echo ','  
	else
		echo ']'
		echo '}'
	fi
 # dump all files concatenated files with the ending in the tmpFile.geojson
done >> tmpFile.geojson
# concatenate the header and the tmpFile to the final file with the name topologies_AgeStart_AgeEndMa.geojson
cat tmpHeader.txt tmpFile.geojson > topologies_${timeStart%.*}-${timeEnd%.*}Ma.geojson
# remove the temporary files
rm tmpFile.geojson tmpHeader.txt 
echo 'Finished'
