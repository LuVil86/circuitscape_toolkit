

import numpy as np
import scipy.ndimage as ndi
import rasterio as rio
from datetime import datetime
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed


def cleanAndFilter(inputRaster,toRemove, toMerge,minPatchSize) :  
    matrix=inputRaster.read(1)
    toRemove=[*toRemove]
    toCross=[*toRemove,*toMerge]


    maskRemove = np.array([[elem in toRemove for elem in row] for row in matrix]) 
    
    
    maskBoth = np.array([[elem in toCross for elem in row] for row in matrix]) 
    maskBoth = np.where(maskBoth==True, 1,0)
    ### le masque de l'érosion donne les endroits où l'algorithme doit opérer (il évitera les autres)
    ## -> on demande de faire l'érosion que sur les chemins
    
    #finMask = ndi.binary_erosion(maskBoth,iterations=5,mask=maskRoad, brute_force=False)
    erosion = ndi.binary_erosion(maskBoth, iterations=5, mask=maskRemove)
    ### le problème c'est que ça érode aussi les chemins dans les patches de forêts QUI TOUCHENT D'AUTRES CLASSES D'HABITATS QUE LES FORÊTS (car inscrits en "0" dans le maskBoth)
    ##  -> donc il faut trouver un moyen de remplacer ces érosions par des valeurs sans pour autant le faire sur les extérieurs des patches
    labeled_array, num_features = ndi.label(erosion, structure=ndi.generate_binary_structure(2,2)) 
    for i in range(1,num_features):
         if np.sum(np.where(labeled_array==i, 1,0)) <= minPatchSize:
              labeled_array[labeled_array==i] = 0
        
    finArray = np.where(labeled_array != 0, 1, 0)
    return finArray



if __name__=="__main__":
    start = datetime.now()
    ######### input parameters ########
    toRemove= [3,24]
    toMerge=[15,13,16,17]
    minPatchSize = 50
    ### input and output
    inputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_15avril25_clip.tif"
    outputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_15avril25_clip_BINARY_PATCHES.tif"


    with rio.open(inputRaster) as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
        cellSize=inp_transform[0]

        result = cleanAndFilter(inp, toRemove, toMerge,minPatchSize )

    out_meta.update({"driver": "GTiff",
                                                "height": result.shape[0],
                                                "width": result.shape[1],
                                                "dtype":"int16",
                                                "nodata":-999
                                                })

    print("-- trying to write raster... ")
    try:
        with rio.open(outputRaster,"w",**out_meta) as dst:
                dst.write(result,1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
