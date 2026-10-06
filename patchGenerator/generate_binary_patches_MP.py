
import numpy as np
import scipy.ndimage as ndi
import rasterio as rio
from datetime import datetime
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
#assert np.__version__>=1.24


def cleanAndFilter(inputRaster,toRemove, toMerge,minPatchSize,x1,x2,y1,y2) :  
    start = datetime.now()
    with rio.open(inputRaster) as src:
        matrix=src.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
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
       # print(f" cleanAndFilter executed in {datetime.now() - start} seconds ")
        return finArray.astype(np.int16)


if __name__=="__main__":
    start = datetime.now()
    ######### input parameters ########
    nbProcessors=16
    toRemove= [3,24]
    toMerge=[15,13,16,17]
    minPatchSize = 12000
    tileSize=(2000,2000)

    ### input and output
    inputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_15avril25_rAoi.tif"
    outputRaster="/home/luvil/test_cleanPathway_fillHoles/results_BINARY_PATCHES_30ha.tif"


    finResults={}
    with rio.open(inputRaster) as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
        cellSize=inp_transform[0]
        nbSplitY=int(np.ceil(nrow/tileSize[0]))
        nbSplitX=int(np.ceil(ncol/tileSize[1]))
           
        print("#"*10)
        print("input raster features : ")
        print(f"{nrow} row X {ncol} columns ")
        print(f"cell size : {cellSize}")
        print(f"'nodata' code : {inp.nodata}")

        print("number of tiles :")
        print(f"{nbSplitY} tiles per row X {nbSplitX} tiles per column")
        if nbSplitY == 1 and nbSplitX ==1 :
            print("!! ERROR only one tile would be generated with the parameters you specified : use the non-multiprocess script instead")
            raise ValueError

        
        splitCoordsXStart=[ i[0] for i in np.array_split(np.arange(ncol),nbSplitX)]  ## number of splits in the x coordinates (so this is the "vertical cuts")
        splitCoordsXStop=[ i[-1] for i in np.array_split(np.arange(ncol),nbSplitX)]
        splitCoordsYStart=[ i[0] for i in np.array_split(np.arange(nrow),nbSplitY)]  ## number of splits in the x coordinates (so this is the "horizontal cuts")
        splitCoordsYStop=[ i[-1] for i in np.array_split(np.arange(nrow),nbSplitY)]
        tileIndex=[ (i,j) for i in range(len(splitCoordsYStart)) for j in range(len(splitCoordsXStart)) ]
    inp.close()


    with ProcessPoolExecutor(max_workers=nbProcessors) as executor:
        poolDF={}


        for k in tileIndex:        
                poolDF[executor.submit(cleanAndFilter,inputRaster,
                                       toRemove,
                                       toMerge,
                                       minPatchSize,
                                       splitCoordsXStart[k[1]],
                                       splitCoordsXStop[k[1]],
                                       splitCoordsYStart[k[0]],
                                       splitCoordsYStop[k[0]],
                                        )]=k
                                       
            
        for future in as_completed(poolDF):
            
            try:
               
               finResults[poolDF[future]]=future.result()
               print(f"tile {poolDF[future]} done")
            except Exception as exc:
                print('%r generated an exception: %s' % (poolDF[future], exc))


######### ***GATHERING RESULTS **** ########

    tmpRow=[]
    for i in range(len(splitCoordsYStart)):
        tmpRow.append(np.hstack(tuple([finResults[(i,j)] for j in range(len(splitCoordsXStart))])))
        
        print(f"column {i} stacking done")
    finResults=None






  #  for i in tmpRow:
   #     print(f"{i.shape[0]} rows  : {i.shape[1]} columns")

    finalMat=np.vstack(tuple([i for i in tmpRow])).astype(np.int16)
    print("row stacking done")
    print("  ##############  ARRAY DONE  ###########")
    print(f" -- final array size : {finalMat.shape[0]} rows X {finalMat.shape[1]} columns")
    print (f" total elapsed time : {datetime.now()-start}")
    
    print("  ##############  ARRAY DONE  ###########")
    print (f" total elapsed time : {datetime.now()-start}")
    out_meta.update({"driver": "GTiff",
                                                "height": finalMat.shape[0],
                                                "width": finalMat.shape[1],
                                                "dtype":"int16",
                                                "nodata":-999
                                                })

    print("-- trying to write raster... ")
    try:
        with rio.open(outputRaster,"w",**out_meta) as dst:
                dst.write(finalMat.astype(np.int16),1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
