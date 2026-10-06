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
import os
#assert np.__version__>=1.24

def computeRSF(enviCoeff, habRasterPath, habCoefs, x1,x2,y1,y2,k):
    
    ##### assign coefficients to habitat classes , stock as first raster ########
    
    with rio.open(habRasterPath) as habPath:
        habVar = habPath.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
        
        habTerm = np.zeros(habVar.shape)
        ###entry coefficients of landuses
        ##create a dictionary for corrrespondances from a csv
        data_dict = pd.Series(habCoefs.coef.values, index=habCoefs.value).to_dict()
        #Calculate the habitat terms b*landuse
        for key, value in data_dict.items():
            habTerm[habVar==key] = value
        TermList = []
        TermList.append(habTerm)

    ###### multiply each continuous variable raster with its respective coefficient  --> add them to TermList

    for covar in enviCoeff:
        with rio.open(covar) as src:
            arr = src.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
            TermList.append((arr/1000)*enviCoeff[covar])
     
    #### stack everything, sum , exponentiate and return #######  
    
    TermStack=np.stack(TermList, axis=0)
    sumTerms = np.sum(TermStack, axis=0)
    expsumTerms = np.exp(sumTerms, dtype=np.float32)     
    return expsumTerms


if __name__=="__main__":
    start = datetime.now()

    ################################ *** input parameters *** ###########################################
        
    habRasterPath='/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Habitats/Hermine/HabitatMap_v1_50_hermine_260106_Nibble2.tif'

    chosenPA = "PA3"
    outFile = '/media/loreto/Grande/ie-ofev-24-25/sdm_hermine/wx_GLM13_AGGR_BCKG_redone.tif'

    ### glm15 = only densities + distEAUX+distPath+landUse
    ### glm16 = only densities + distEAUX+distPath
    GLMcoeff = pd.read_csv('/media/loreto/T7_ROUGE/lolo/sdm_opportunistic_data/RSF_with_Reass_2_output/glm13_AGGR_BCKG_50m_nbPA_100000_nb_Draw_5.csv')
    GLMcoeff["new_short_class"] = GLMcoeff["term"].str.replace("landUse", "")
    GLMcoeff.rename(columns={"Estimate" : "coef"}, inplace=True)

    NASpath= "/media/loreto/NAS_DEVELOPPEMENT"

    covarRasterPath=pd.read_csv('/media/loreto/Grande/ie-ofev-24-25/sdm_hermine/sdm_hermine_covariable_path_list.csv')
    covarRasterPath["rasterPath"] = NASpath + covarRasterPath["rasterPath"].astype(str)
 

    ######### * COEFFICIENTS LAND USE * #########
    
    #### coefs du glm global (landUse + variables continues)
    
    correspHabitat = pd.read_csv(os.path.join(NASpath,'IE_OFEV/Habitats/Hermine/unique_habitat_hermine_for_resistances_06janv26_with_expert_ranking.csv'))
    habCoefs = pd.merge( correspHabitat,GLMcoeff, on="new_short_class")
    habCoefs = habCoefs.loc[habCoefs["PA-set"] == chosenPA]

    ### avis d'expert

    '''
    habCoefs = pd.read_csv('/media/loreto/Grande/ie-ofev-24-25/sdm_hermine/res_expert/unique_habitat_hermine_for_resistances_06janv26_with_expert_ranking.csv')
   # habCoefs["log_coef"] = np.log(habCoefs['exp_ranking'])
   # habCoefs.rename(columns={'log_coef':'coef'},inplace=True)

    #### "exp_ranking3" is for resistances : "exp_ranking" is for modelling w(x)
    habCoefs.rename(columns={'exp_ranking':'coef'},inplace=True)
    '''

    ### coeff nuls (pour générer un modèle avec uniquement les covariables continues)
    '''
    habCoefs = pd.read_csv('/media/loreto/Grande/ie-ofev-24-25/sdm_hermine/res_expert/unique_habitat_hermine_for_resistances_06janv26_with_expert_ranking.csv')
    habCoefs["coef"] = 0
    '''
    ### coeff GLM effectué avec landUse only
    #GLMcoeffHab = pd.read_csv('/media/luvil/T7_ROUGE/lolo/sdm_opportunistic_data/RSF_with_Reass_2_output/glm12_AGGR_BCKG_50m_nbPA_100000_nb_Draw_5.csv')
    #GLMcoeffHab["new_short_class"] = GLMcoeffHab["term"].str.replace("landUse", "")
    #GLMcoeff.rename(columns={"Estimate" : "coef"}, inplace=True)
    #correspHabitat = pd.read_csv('/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/Habitats/Hermine/unique_habitat_hermine_for_resistances_02janv26.csv')
    #habCoefs = pd.merge( correspHabitat,GLMcoeff, on="new_short_class")
    #habCoefs = habCoefs.loc[habCoefs["PA-set"] == chosenPA]

    #####################################################################################



    

    structCoeff =pd.merge(GLMcoeff,covarRasterPath, left_on="new_short_class", right_on="rasterName")
    structCoeff = structCoeff.loc[structCoeff["PA-set"] == chosenPA]
    structCoeff = structCoeff.loc[structCoeff["Pr(>|z|)"] < 0.05]
    enviCoeff =pd.Series(structCoeff.coef.values, index=structCoeff.rasterPath).to_dict()
    
    print("****** CONTINUOUS VARIABLE COEFFICIENTS ****** ")
    print(structCoeff[["rasterName","coef"]])
    print("****** LAND USE COEFFICIENTS ****** ")
    print(habCoefs[["value", "new_short_class", "coef"]])
    
                       
    tileSize=(6000,6000)
    nbProcessors = 10
        
        
        
    
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
        print(f"dtype : {inp.dtypes}")
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
            tmpRow.append(np.hstack(tuple([finResults[(i,j)] for j in range(len(splitCoordsXStart))])))
            
            print(f"column {i} stacking done")
            
    ### remove MP results #####
    finResults=None



    finalMat=np.vstack(tuple([i for i in tmpRow]))
    print("row stacking done")
    print("  ##############  ARRAY DONE  ###########")
    print(f" -- final array size : {finalMat.shape[0]} rows X {finalMat.shape[1]} columns")
    print (f" total elapsed time : {datetime.now()-start}")
    print(f"smaller value of minExp : {min(minExp)}, bigger value of maxExp : {max(maxExp)} ")
    
    out_meta.update({"driver": "GTiff",
                        "height": finalMat.shape[0],
                        "width": finalMat.shape[1],
                        "dtype":"float32",
                        "nodata":np.nan
                        })

    print("-- writing raster... ")
    try:
        with rio.open(outFile,"w",**out_meta) as dst:
                dst.write(finalMat,1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
        print(" >>>> SCRIPT COMPLETED")
