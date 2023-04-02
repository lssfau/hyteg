---------------------------

Input files for calculating the surface velocities were obtained or compiled
from Müller R.D., Seton, M., Zahirovic, S., Williams, S.E., Matthews, K.J.,
Wright, N.M., Shephard, G.E., Maloney, K.T., Barnett-Moore, N., Hosseinpour,
M., Bower, D.J., Cannon, J., 2016. Ocean basin evolution and global-scale
plate reorganization events since Pangea breakup, Annual Review of Earth and
Planetary Sciences, Vol 44, 107-138. DOI: 10.1146/annurev-earth-060115-012211.

Files used: 
Global_EarthByte_230-0Ma_GK07_AREPS.rot
Global_EarthByte_230-0Ma_GK07_AREPS_PlateBoundaries.gpml

---------------------------
Plates input files
---------------------------
The rotational is the exact same file - Global_EarthByte_230-0Ma_GK07_AREPS.rot
The plate boundaries file (topologies.geojson) is a pre-processed file. It
contains the plate configuration for each desired time 0-100Ma with a time
step of 0.1Ma. This file is a merge of all the geojson files extracted from
GPlates using the rotational and the plate boundaries gpml from the global
model by Müller et al. (2016). Each collection has a name tag with the age
formatted to 4 digits, e.g. 0.0000Ma.

To extract such files from GPlates:
Reconstruction $>$ Extraction $>$ Add Export $>$ 

1. Resolved Topologies (General) 
2. GeoJSON (*.geojson) format
3. Tick in "export into a single file" and "export resolved topological
   polygons". Unclick the rest.  
4. File name as "topology\%P\_\%0.4fMa".

GPlates extracts the plate boundaries in a separate file for each age required.
To use it in the HyTeG framework, the files were groupped into one file using
a bash script, following the json format. 

---------------------------
Contact details: 
Berta Vilacís - bvilacis@geophysik.uni-muenchen.de
---------------------------
