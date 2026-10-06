
import numpy as np
import scipy.ndimage as ndi
import rasterio as rio
from datetime import datetime
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
#assert np.__version__>=1.24

def clean_pathway(inputRaster,toRemove, milieux,minPatchSize,nNeighbor,x1,x2,y1,y2,padding,nrow,ncol,k) : 

    nrowOrig=y2-y1
    ncolOrig=x2-x1
    toRemove=[*toRemove]
    milieux=[*milieux]
#  print(f"y1 = {y1}, y2 = {y2}, x1 = {x1}, x2={x2}")
#    print(f"nrowOrig = {nrowOrig} , ncolOrig = {ncolOrig}")
    if (x1-padding<0 or y1-padding<0 or x2+padding>ncol or y2+padding>nrow):  
        print(f"the tile {k} is a border tile : no computation required")
        return np.zeros((nrowOrig,ncolOrig), dtype=np.uint8)
    else:
        with rio.open(inputRaster, 'r') as rasterBuffer:  
            add_mat = nNeighbor  
            matrix=rasterBuffer.read(1,window=Window.from_slices((y1-padding, y2+padding), (x1-padding, x2+padding)),out_dtype=np.uint8)

            inflatedZero= np.pad(np.zeros(matrix.shape) ,((nNeighbor,nNeighbor), (nNeighbor,nNeighbor)), mode="constant", constant_values=0)
            inflatedMat = np.pad(np.zeros(matrix.shape) ,((nNeighbor,nNeighbor), (nNeighbor,nNeighbor)), mode="constant", constant_values=0)
              

            #### generate inflated matrix : split into the two elements                                           
            inflatedMat[add_mat:-add_mat,add_mat:-add_mat] = matrix                 
            maskRemove = np.array([[elem in toRemove for elem in row] for row in inflatedMat],dtype=np.uint8)                                          
            maskMilieux = np.array([[elem in milieux for elem in row] for row in inflatedMat],dtype=np.uint8)  


            for i in range(add_mat,inflatedMat.shape[0]-add_mat) :                          
                for j in range(add_mat,inflatedMat.shape[1]-add_mat):                       
                           
                    if maskRemove[i,j] and np.any(maskMilieux[i-add_mat:i+add_mat+1, j-add_mat:j+add_mat+1]):                                  
                        inflatedZero[i,j] = 1                                              
            maskMilieux[inflatedZero==1] = True
            tmp=maskMilieux[add_mat:-add_mat,add_mat:-add_mat]
            erosion = np.where(tmp, 1, 0)
            
        labeled_array, num_features = ndi.label(erosion, structure=ndi.generate_binary_structure(2,2)) 

        #### method to filter patches with "square" delineation : the patches are filtered by their extent size. This is faster, but linear patches size
        #### would be heavily overestimated and might pass the filter.
        for f in ndi.find_objects(labeled_array):
           dim_f =  (f[0].stop-f[0].start)*(f[1].stop-f[1].start)
           if dim_f<= minPatchSize:
               labeled_array[f] = 0
        #### method to filter patches with actual number of pixels in the patch. this is the most accurate way in terms of surface, but is very slow
  #      for i in range(1,(num_features+1)):
  #          if np.sum(np.where(labeled_array==i, 1,0)) <= minPatchSize:
  #              labeled_array[labeled_array==i] = 0


        finArray = np.where(labeled_array[padding:padding+nrowOrig,padding:padding+ncolOrig] != 0, 1, 0)
       # print(f" Nrow : {finArray.shape[0]} x Ncol : {finArray.shape[1]}")
        return finArray.astype(np.uint8)




if __name__=="__main__":
    start = datetime.now()
    ######### input parameters ########
    nbProcessors=16
    toRemove= [3,21,24]
    milieux=[15,13,16,17]
    '''
    inputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_15avril25_clip.tif"
    minPatchSize = 100
    tileSize=(100,100)
    '''
    percentBuffer=200
    nNeighbor=1
    ### input and output
    
    inputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_forZN_rAoi.tif"
    minPatchSize = 16000
    tileSize=(1000,1000)
    
    outputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_forZN_cleanPathway_BINARY_PATCHES_w_secRoad.tif"


    finResults={}
    with rio.open(inputRaster) as inp:
        expFactor=percentBuffer/100
        padding=int(tileSize[0]*expFactor) ### the number of cells to add on the contours
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
                poolDF[executor.submit(clean_pathway,inputRaster,
                                       toRemove,
                                       milieux,
                                       minPatchSize,
                                       nNeighbor,
                                       splitCoordsXStart[k[1]],
                                       splitCoordsXStop[k[1]]+1,
                                       splitCoordsYStart[k[0]],
                                       splitCoordsYStop[k[0]]+1,
                                        padding,
                                        nrow,
                                        ncol
                                        ,k)]=k
                                                
            
        for future in as_completed(poolDF):
            
            try:
               
               finResults[poolDF[future]]=future.result()
               print(f"tile {poolDF[future]} done")
            except Exception as exc:
                print('%r generated an exception: %s' % (poolDF[future], exc))


######### ***GATHERING RESULTS **** ########

    tmpRow=[]
    for i in range(len(splitCoordsYStart)):
        tmpRow.append(np.hstack(tuple([finResults[(i,j)] for j in range(len(splitCoordsXStart))]),dtype=np.uint8))
        
        print(f"column {i} stacking done")
    finResults=None






  #  for i in tmpRow:
   #     print(f"{i.shape[0]} rows  : {i.shape[1]} columns")

    finalMat=np.vstack(tuple([i for i in tmpRow])).astype(np.uint8)
    print("row stacking done")
    print("  ##############  ARRAY DONE  ###########")
    print(f" -- final array size : {finalMat.shape[0]} rows X {finalMat.shape[1]} columns")
    print (f" total elapsed time : {datetime.now()-start}")
    
    print("  ##############  ARRAY DONE  ###########")
    print (f" total elapsed time : {datetime.now()-start}")
    out_meta.update({"driver": "GTiff",
                                                "height": finalMat.shape[0],
                                                "width": finalMat.shape[1],
                                                "dtype":np.uint8,
                                                "nodata":0
                                                })

    print("-- trying to write raster... ")
    try:
        with rio.open(outputRaster,"w",**out_meta) as dst:
                dst.write(finalMat,1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
