
import numpy as np
import pandas as pd
from datetime import datetime

import geopandas as gp
import sys
from rasterio.features import shapes
from shapely.geometry import shape, Polygon
from rasterio.windows import Window

import rasterio as rio
import json as jsn
from scipy import ndimage as ndi
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,40).__str__()
import cv2 as cv
import matplotlib.pyplot as plt
import geopandas as gpd


def polygonize(fileReader,numPatch):
    time=datetime.now()
    shapeObject=ndi.find_objects(fileReader.read(1)==numPatch )
    win=Window.from_slices(shapeObject[0][0], shapeObject[0][1])
    subRast=fileReader.read(1,window=win, out_dtype=np.int32)
    results = (
    {'raster_val':v, 'geometry': s}
    for i, (s, v) 
    in enumerate(
    shapes(source=subRast,transform=fileReader.transform, connectivity=8)))
    test=list(filter(lambda x:x["raster_val"]==numPatch,results))
    fin=gpd.GeoSeries(shape(test[0]["geometry"]))
    print(f"total elapsed time : {datetime.now()-time}")
    return fin
if __name__=='__main__':
    numPatch=63188
    with rio.open("/home/lucas/patches_19_20_21_10m.tif") as src:
        maxPatchNum=int(src.statistics(1).max)
        print(maxPatchNum)
        patchList=polygonize(src, numPatch=numPatch)
   
    