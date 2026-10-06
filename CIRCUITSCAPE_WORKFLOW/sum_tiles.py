#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 18:35:16 2023

@author: loreto
"""

### This should calculate the sum of a lot of rasters as 
##long as they are all the same resolution, extent, CRS, 
##etc. I put in an assert statement to double check that.


import rasterio
import os

# ******** INPUT PARAMETERS ******** #
cutTileFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_jura/cut_tuile_20"

outFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_jura/summed_tuile_20"
# *********************************** #


if __name__=='__main__':
    print("******** SUM EAST-WEST AND NORTH-SOUTH TILES FROM CIRCUITSCAPE OUTPUTS *******\n")
    print(f"folder with tiles : {cutTileFolder}")
    print(f"output folder : {outFolder}")
    print("*"*100)
    tileList=[]
    
    if not os.path.exists(outFolder):
	    print("the specified output folder does not exists :: creating it")
	    os.makedirs(outFolder)
		
    for file in os.listdir(cutTileFolder):
        if file.endswith("ns_cum_curmap.tif"):
            tileList.append(int(file.split("_")[2]))
            
    for tileNumber in tileList:
        matchList=[]    
        dimList = []
        for file in os.listdir(cutTileFolder):
            if file.find(f"_{tileNumber}_")!=-1:
                with rasterio.open(os.path.join(cutTileFolder,file)) as src:    
                    tmp = src.read(1)
                    dimList.append(tmp.shape)
                    matchList.append(tmp)
                    result_profile=src.profile
        if (len(matchList) == 2) and (dimList[0] == dimList[1]):
            try:
                with rasterio.open(os.path.join(outFolder, f't_{tileNumber}.tif'), 'w', **result_profile) as dst:
                    rasterFinal=matchList[0]+matchList[1]
                    dst.write(rasterFinal,1)            
                    print(f"tile {tileNumber} done")
            except Exception as e:
                print(f"WARNING : there was an error writing summed tile N°{tileNumber} : skipping this tile...")
                print(e)
                continue
        else:
            print(f"""WARNING : the number of CircuitScape tile results for the tile {tileNumber} is not 2 or the two 
            tiles are not the same size : skipping the tile""")
            continue
    print(" >>>> SCRIPT COMPLETED")
        
