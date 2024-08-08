import os
import json
import numpy as np
import traceback


rotFileIn = '../Datasets/Chen2023/EB_410to0Ma_GK07_Matthews_2016_v2.rot'
rotFileIn = '../Datasets/Chen2023/AR_410to0Ma.rot'
rotFileIn = '../Datasets/Muller_etal_2022_SE_1Ga_Opt_PlateMotionModel/Muller_etal_2022_SE_1Ga_Opt_PlateMotionModel/optimisation/no_net_rotation_model.rot'
rotFileIn = '../Datasets/Gibbons_etal_EarthByte_PlateModel_GPlates/Global_EarthByte_TPW_CK95G94_Rigid_Gibbons.rot'


# we are going to read each line and append longitude and latitudes
plateID = []
time = []
lat = []
lon = []
angle = []
conjPlateID = []
comment = []
#reading rotational
with open(rotFileIn, mode='r') as my_id:
    for my_line in my_id:
        try:
            splt_my_line = my_line.split()
            plateID.append(int(splt_my_line[0]))
            time.append(float(splt_my_line[1]))
            lat.append(float(splt_my_line[2]))
            lon.append(float(splt_my_line[3]))
            angle.append(float(splt_my_line[4]))
            conjPlateID.append(int(splt_my_line[5]))
            comment.append(' '.join(splt_my_line[6:len(splt_my_line)+1]))
            
        except Exception as e:
            print(my_line)
            print('----------------------')
            traceback.print_exc()
            print()
            print('----------------------')
            raise e

rotData = np.array((plateID, time, lat, lon, angle, conjPlateID, comment), dtype = 'object').transpose()
rotData = rotData[rotData[:,1].argsort()] # sort by time
rotData = rotData[rotData[:,0].argsort(kind='mergesort')] # sort by PlateID
fmt = "%03d %0.1f %0.4f %0.4f %0.4f %03d %-40s"
np.savetxt(rotFileIn[:-4]+"-ordered.rot", rotData, fmt = fmt, delimiter = "  ")
print('File written in ::'+rotFileIn[:-4]+"-ordered.rot")

# Swap function
def swapPositions(arr, pos1, pos2):
    arr[pos1], arr[pos2] = arr[pos2].copy(), arr[pos1].copy()
    return arr

IndexToBeDeleted = []
rotDataFiltered = rotData.copy()
for idx, line in enumerate(rotData): 
    if  idx < len(rotData)-1:
        line2 = rotData[idx+1]
        if(len(line)==len(line2) and len(line)-1==sum([1 for i,j in zip(line[0:-1],line2[0:-1]) if i==j])):
            IndexToBeDeleted.append(idx+1)
        elif line[0] == line2[0] and line[1] == line2[1]: 
            if line[5] != rotData[idx-1][5] and line2[5] != rotData[idx+2][5]:
                swapPositions(rotDataFiltered, idx, idx+1)

rotDataFiltered= np.delete(rotDataFiltered, IndexToBeDeleted, axis = 0)
# Consider special case when one has four lines (two and two) with the same age and plate ID. 
# Sort them accordingly. 
for idx, line in enumerate(rotDataFiltered): 
    if  idx < len(rotDataFiltered)-3:
        line2 = rotDataFiltered[idx+1]
        line3 = rotDataFiltered[idx+2]
        line4 = rotDataFiltered[idx+3]

        if line[0] == line2[0] and line[1] == line2[1] and line3[0] == line4[0] and line3[1] == line4[1]: 
            if line[5] != rotDataFiltered[idx-1][5] and line2[5] == line3[5] and line4[5] != rotDataFiltered[idx+4][5]:
                swapPositions(rotDataFiltered, idx, idx+1)
                swapPositions(rotDataFiltered, idx+2, idx+3)


print( rotData.shape, rotDataFiltered.shape)

fmt = "%03d %0.1f %0.4f %0.4f %0.4f %03d %-40s"
np.savetxt(rotFileIn[:-4]+"-sorted.rot", rotDataFiltered, fmt = fmt, delimiter = "  ")
print('File written in ::'+rotFileIn[:-4]+"-sorted.rot")