#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 11:33:57 2025

@author: loreto
"""


import numpy as np
import pandas as pd
import rasterio as rio
from datetime import datetime
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
#assert np.__version__>=1.24

def computeRSF(enviCoeff, habRasterPath, habCoefs, x1,x2,y1,y2,k):
    #print(f" running resistance calculation for tile {k}")  
##### assign coefficients to habitat classes ########
    with rio.open(habRasterPath) as habPath:
        habVar = habPath.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
        
        habTerm = np.empty(habVar.shape, dtype=np.dtype('float32'))
        ###entry coefficients of landuses
        ##create a dictionary for corrrespondances from a csv
        data_dict = habCoefs.set_index('value').to_dict()['coef']
        #Calculate the habitat terms b*landuse
        for key, value in data_dict.items():
            habTerm[habVar==key] = value
### stock habitat coefficient raster within termList
    TermList=[]
    TermList.append(habTerm)

###### multiply each continuous variable raster with its respective coefficient 
###### --> add them to TermList    
    for covar in enviCoeff:
        with rio.open(covar) as src:
            arr=src.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
            TermList.append((arr/1000)*enviCoeff[covar])
     
    #### stack everything #######   
            TermStack=np.stack(TermList, axis=0)
            sumTerms = np.sum(TermStack, axis=0, dtype=np.float32)
            expsumTerms = np.exp(sumTerms, dtype=np.float32)
    return expsumTerms

def computeResistance(expRaster, minExp, maxExp):
    hsI=(expRaster-minExp)/(maxExp-minExp)    
    kterm=np.exp(-1*4*hsI, dtype=np.float32)
    res = 100-99*((1-kterm)/(1-np.exp(-4, dtype=np.float32)))
    return res

if __name__=="__main__":
    start = datetime.now()

    ################################ *** input parameters *** ###########################################
        
    habRasterPath='/media/loreto/Grande/ie-ofev-24-25/variables/habitat_cerf_for_resistances/HabitatMap_cerf_07juil25_for_resistances.tif'
    
    habCoefs= pd.read_csv('/media/loreto/Grande/ie-ofev-24-25/ssf_jura/ssf_raster/unique_habitat_cerf_for_resistance_coef_22_model.csv')
    
    enviCoeff={'/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/Jura__Density_Buildings_100_opt2_scaled.tif':0.1569276,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/Jura__Density_Merge_RoadPrimary__200_opt_scaled.tif':0.1541613,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/Jura__Density_Merge_RoadSecondary_50_opt_scaled.tif':-0.3041774,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/Jura__Density_Forest_400_opt2_scaled.tif':0.3176758,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/Jura__Dist_Merge_Bati_16b_scaled.tif':0.1605695,
               #'/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/Jura__Dist_Merge_RoadPrimary_16b_scaled.tif':0.0964482,
               #'/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/Jura__Dist_Merge_Autobahn_16b_scaled.tif':0.0901117,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/Jura__Dist_Merge_RoadSecondary_250701_16b_scaled.tif':-0.3298851,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/Scaled/Jura__Altitude_5m_16b_scaled.tif':0.3157216,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/Scaled/Jura__Exposition_5m_16b_scaled.tif':-0.0571324,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/Scaled/Jura__Slope_5m_8b_scaled.tif':-0.2298444
               }
                     
                       
    tileSize=(6000,6000)
    nbProcessors = 10
        
        
        
    outFile = '/media/loreto/Grande/ie-ofev-24-25/ssf_jura/ssf_raster/output/resistance_cerf_jura_19_sep_25_model20.tif'
    
    ############################################################################################3
        
    finResults={}
    with rio.open(habRasterPath) as inp:
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
    
    
            
        splitCoordsXStart=[ i[0] for i in np.array_split(np.arange(ncol),nbSplitX)] 
        splitCoordsXStop=[ i[-1] for i in np.array_split(np.arange(ncol),nbSplitX)]
        splitCoordsYStart=[ i[0] for i in np.array_split(np.arange(nrow),nbSplitY)]
        splitCoordsYStop=[ i[-1] for i in np.array_split(np.arange(nrow),nbSplitY)]
        tileIndex=[ (i,j) for i in range(len(splitCoordsYStart)) for j in range(len(splitCoordsXStart)) ]
        inp.close()
    
    with ProcessPoolExecutor(max_workers=nbProcessors) as executor:
            poolDF={}
            minExp=[]
            maxExp=[]
            
            for k in tileIndex:        
                    poolDF[executor.submit(computeRSF,
    									   enviCoeff,
                                           habRasterPath,
                                           habCoefs,
                                           splitCoordsXStart[k[1]],
                                           splitCoordsXStop[k[1]],
                                           splitCoordsYStart[k[0]],
                                           splitCoordsYStop[k[0]],
                                           k)]=k
                                           
            for future in as_completed(poolDF):
                try:
                   finResults[poolDF[future]]=future.result()
                   minExp.append(np.min(finResults[poolDF[future]]))
                   maxExp.append(np.max(finResults[poolDF[future]]))
                   print(f"tile {poolDF[future]} done")
                except Exception as exc:
                    print('%r generated an exception: %s' % (poolDF[future], exc))
    
    
    ######### ***GATHERING RESULTS **** ########
        
        
    tmpRow=[]
    smallerExp=min(minExp)
    biggerExp=max(maxExp)
    for i in range(len(splitCoordsYStart)):
            tmpRow.append(np.hstack(tuple([computeResistance(finResults[(i,j)], smallerExp,biggerExp) for j in range(len(splitCoordsXStart))])))
            
            print(f"column {i} stacking done")
            
    ### remove MP results #####
    finResults=None



    for i in tmpRow:
        print(f"{i.shape[0]} rows  : {i.shape[1]} columns")
       
        finalMat=np.vstack(tuple([i for i in tmpRow])).astype(np.int16)
        
        print("row stacking done")
        print("  ##############  ARRAY DONE  ###########")
        print(f" -- final array size : {finalMat.shape[0]} rows X {finalMat.shape[1]} columns")
        print (f" total elapsed time : {datetime.now()-start}")
        print(f"smaller value of minExp : {min(minExp)}, bigger value of maxExp : {max(maxExp)} ")
        
        out_meta.update({"driver": "GTiff",
                         "height": finalMat.shape[0],
                         "width": finalMat.shape[1],
                         "dtype":"int16",
                         "nodata":-9999
                         })
    
    print("-- writing raster... ")
    try:
        with rio.open(outFile,"w",**out_meta) as dst:
                dst.write(finalMat,1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
        print(" >>>> SCRIPT COMPLETED")
