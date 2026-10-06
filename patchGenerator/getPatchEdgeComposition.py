

import numpy as np
import pandas as pd
from datetime import datetime
import pylandstats as pls
import geopandas as gp
import sys
import rasterio as rio
from rasterio.mask import mask
import json as jsn




def getPatchEdgeComposition(patchMapFile, habitatMapFile, unit="pixel"):
    start = datetime.now()
    habitatMap=rio.open(habitatMapFile,'r')
    print("===> habitat map loaded.")

    patches=gp.read_file(patchMapFile) 
    patches.sort_values(by="patches", inplace=True)
    patches=patches.reset_index(drop=True)

          #read the shp with the forest patches
    print(" ===> patches file loaded")

    def getFeatures(gdf):
        return [jsn.loads(gdf.to_json())['features'][0]['geometry']]
    
    DFinit=pd.DataFrame(index=np.unique(habitatMap.read(1)))

    tot=len(patches)
    finDict={}
    for i in range(tot):

        coords_original=gp.GeoDataFrame({"geometry":patches.loc[i].geometry}, index=[0])
        buffed=patches.iloc[i].geometry.buffer(distance=10, resolution=1, join_style=2)     #create buffer around the patch (10m = size of the pixel)
        coords=gp.GeoDataFrame({"geometry":buffed}, index=[0])
        maskedHabitat, transformed =mask(habitatMap,getFeatures(coords), crop=True, nodata=0) 
       # maskedHabitat=np.squeeze(maskedHabitat)
       # ls = pls.Landscape(maskedHabitat, res=(10,10)) 
                #clip the habitat map with the buffer
        out_meta=habitatMap.meta
        out_meta.update({"driver":"GTiff",                                                  #need to write the maskedhabitat to re-clip the patch on it
                        "height":maskedHabitat.shape[1],
                        "width":maskedHabitat.shape[2],
                        "transform":transformed,
                        "nodata": 0})                                                    #define no data different from -128 to discriminate between real no data and artefact no data created by the mask
        with rio.open("tmp_mask", "w", **out_meta) as final:
            final.write(maskedHabitat)
        with rio.open("tmp_mask") as maskedHabitat:
            mask2, transf2=mask(maskedHabitat,getFeatures(coords_original),invert=True, nodata=0) 
        mask2=np.squeeze(mask2)     #squeeze remove one dimension
        ls = pls.Landscape(mask2, res=(10,10))         #ls is an object with Landscape type  res=(mask2.shape[0],mask2.shape[1])
        if unit == "pixel":
            class_area=ls.compute_class_metrics_df(metrics=["total_area"], metrics_kws={"total_area":{"hectares":False}})  #metrics_kws={"total_area":{"hectares":False}
        elif unit =="hectare":
            class_area=ls.compute_class_metrics_df(metrics=["total_area"], metrics_kws={"total_area":{"hectares":True}})  #metrics_kws={"total_area":{"hectares":False}
        else:
            print(" the unit you specified is not valid : choose between 'hectare' or 'pixel' ")
            sys.exit(1)
        class_area.rename(columns={'total_area':str(i+1)},inplace=True)   
        finDict[i]=pd.concat([DFinit, class_area], axis=1)

        print(f"Patch {i+1}/{tot} done", end="\r")

    finDF=pd.concat([v for v in finDict.values()], axis=1)
    print(" ======= DONE =========")
    print("Time for the all process: ", datetime.now()-start)
    return finDF.transpose()

if __name__=='__main__':
    patchDF=getPatchEdgeComposition(patchMapFile="/home/lucas/patches_20_21_10m_polygons.shp",
                              habitatMapFile="/home/lucas/PARATAGE CERF - Jacinthe/Habitat_Raster_CLEAN_10m/_120_Hab_Rast_10m.tif",
                              unit="pixel")
    patchDF.to_csv("patch_statistics_test_pixel.csv")