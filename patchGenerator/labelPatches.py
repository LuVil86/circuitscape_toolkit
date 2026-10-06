# -*- coding: utf-8 -*-

import numpy as np
import os
from datetime import datetime
import rasterio as rio
from scipy import ndimage as ndi
    
   

if __name__=="__main__":
    print(" ********* LABEL PATCHES ********")
    ##### parameters
    inputRasterPath="/media/loreto/Grande/ie-ofev-24-25/variables/habitat_cerf_forZN/ZN_cerf_rAoi_BINARY_PATCHES_trail2_wo_clearing_secRoads.tif"
    basename=inputRasterPath.split(".")[0]
    #outputRasterPath="/media/luvil/DATA/HabitatMap_cerf_forZN_cleanPathway_BINARY_PATCHES_w_secRoad_1000_buff100_LABELED.tif"
    outputRasterPath=basename+'_LABELED.tif'



    with rio.open(inputRasterPath) as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
        cellSize=inp_transform[0]

        r = inp.read(1)
        print("-- raster reading done")
        labeled_array, num_features = ndi.label(r, structure=ndi.generate_binary_structure(2,2))
        print(f" {num_features} patches found and labeled")
        out_meta.update({"driver": "GTiff",
                                                    "height": labeled_array.shape[0],
                                                    "width": labeled_array.shape[1],
                                                    "dtype":"int16",
                                                    "nodata":0
                                                    })

        print("-- trying to write raster... ")
        try:
            with rio.Env(CHECK_DISK_FREE_SPACE="NO"):
                with rio.open(outputRasterPath,"w",**out_meta) as dst:
                        dst.write(labeled_array.astype(np.int16),1)   
        except Exception as exc:
            print("writing raster generated an exception : ", exc)
