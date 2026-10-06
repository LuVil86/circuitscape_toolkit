
import numpy as np
import os
import scipy.ndimage as ndi
import rasterio as rio
from rasterio.merge import merge
from datetime import datetime
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
import glob
#assert np.__version__>=1.24


def cleanAndFilter(inputRaster,toRemove, toMerge,minPatchSize,minPatchMethod,x1,x2,y1,y2,padding,nrow,ncol,k, outputTileFolder) :  
    start = datetime.now()
    
    nrowOrig=y2-y1
    ncolOrig=x2-x1
#  print(f"y1 = {y1}, y2 = {y2}, x1 = {x1}, x2={x2}")
#    print(f"nrowOrig = {nrowOrig} , ncolOrig = {ncolOrig}")
    if (x1-padding<0 or y1-padding<0 or x2+padding>ncol or y2+padding>nrow):  ### this is for the boundary tiles
            print(f"the tile {k} is a border tile : no computation required")
            return("BORDER")
    else:
        with rio.open(inputRaster) as src:
            
            in_meta = src.meta.copy()
            win=Window.from_slices((y1-padding, y2+padding), (x1-padding, x2+padding))
            matrix=src.read(1,window=win)


            win_transform = src.window_transform(win)

            in_meta.update({ "height":matrix.shape[0], 
                            "width":matrix.shape[1],
                            "nodata":0,
                            "transform":win_transform})

                
            toRemove=[*toRemove]
            toCross=[*toRemove,*toMerge]
            maskRemove = np.array([[elem in toRemove for elem in row] for row in matrix], dtype=np.uint8) 
            maskBoth = np.array([[elem in toCross for elem in row] for row in matrix], dtype=np.uint8) 
            maskBoth = np.where(maskBoth==True, 1,0)
            ### le masque de l'érosion donne les endroits où l'algorithme doit opérer (il évitera les autres)
            ## -> on demande de faire l'érosion que sur les chemins
            
            #finMask = ndi.binary_erosion(maskBoth,iterations=5,mask=maskRoad, brute_force=False)
            erosion = ndi.binary_erosion(maskBoth, iterations=5, mask=maskRemove)
            ### le problème c'est que ça érode aussi les chemins dans les patches de forêts QUI TOUCHENT D'AUTRES CLASSES D'HABITATS QUE LES FORÊTS (car inscrits en "0" dans le maskBoth)
            ##  -> donc il faut trouver un moyen de remplacer ces érosions par des valeurs sans pour autant le faire sur les extérieurs des patches
            if minPatchSize != 0:
                labeled_array, num_features = ndi.label(erosion, structure=ndi.generate_binary_structure(2,2)) 
                if minPatchMethod == "envelope":
                        for f in ndi.find_objects(labeled_array):
                            dim_f =  (f[0].stop-f[0].start)*(f[1].stop-f[1].start)
                            if dim_f<= minPatchSize:
                                labeled_array[f] = 0
                elif minPatchMethod == "pixel_count":
                    for i in range(1,(num_features+1)):
                        if np.sum(np.where(labeled_array==i, 1,0)) <= minPatchSize:
                            labeled_array[labeled_array==i] = 0
                else:
                    print("method to estimate patch size does not exist: aborting script")
                    raise ValueError
                finArray = np.where(labeled_array != 0, 1, 0)
                with rio.open(os.path.join(outputTileFolder,"Tile_"+str(k)+".tif"), "w",**in_meta) as dest:
                     dest.write(finArray,1)
                     return("OK")
        # print(f" cleanAndFilter executed in {datetime.now() - start} seconds "
            else:
                with rio.open(os.path.join(outputTileFolder,"Tile_"+str(k)+".tif"), "w",**in_meta) as dest:
                     dest.write(erosion,1)
                     return("OK")


if __name__=="__main__":
    start = datetime.now()
    ######### input parameters ########

    nbProcessors=16
    #### habitat to dissolve in the patches
    toRemove= [3,21,24]
    #### habitat that are merged to consider a nodal zone
    toMerge=[13,15,16,17]
    ### minimum nodal zone size in the tile + buffer (in pixel size)
    minPatchSize = 4000    #### 12000 *25sq.m = 300'000sq.m = 30ha :: but beware that it is dependent on tile size therefore big ZN might be "cut and discarded"
    minPatchMethod = "pixel_count"
    ### tile size
    tileSize=(2000,2000)
    ### percent of tile size used to buffer (8 square )
    percentBuffer=100
    
    ### input and output
    inputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_forZN_rAoi.tif"
    outputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_forZN_rAoi_BINARY_PATCHES_2000_buff100_erosion_30Ha_cleaned_secRoad.tif"
    outputTileFolder = "/home/luvil/test_cleanPathway_fillHoles/ZN_tiles/"


    finResults={}
    with rio.open(inputRaster) as inp:
        expFactor=percentBuffer/100
        padding=int(tileSize[0]*expFactor) ### the number of cells to add on the contours
        out_meta=inp.meta.copy()
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
                                       minPatchMethod,
                                       splitCoordsXStart[k[1]],
                                       splitCoordsXStop[k[1]]+1,
                                       splitCoordsYStart[k[0]],
                                       splitCoordsYStop[k[0]]+1,
                                        padding,
                                        nrow,
                                        ncol
                                        ,k, outputTileFolder)]=k
                                                
            
        for future in as_completed(poolDF):
            
            try:
            
                finResults[poolDF[future]]=future.result()
                print(f"tile {poolDF[future]} done")

            except Exception as exc:
                print('%r generated an exception: %s' % (poolDF[future], exc))


    m, out_transform = merge(glob.glob(outputTileFolder+"*.tif"), method="max")
    finArray = np.squeeze(m, axis=0)
    print(finArray.shape)
    out_meta.update({"height":finArray.shape[0],
                      "width":finArray.shape[1],
                      "transform" : out_transform})
    
    with rio.open(outputRaster, "w", **out_meta) as dest:
        dest.write(finArray, 1)