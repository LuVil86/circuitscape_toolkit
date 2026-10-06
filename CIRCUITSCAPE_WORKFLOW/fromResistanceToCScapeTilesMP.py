import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import rasterio as rio
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor,ThreadPoolExecutor, as_completed
import sys
from datetime import datetime
from pathlib import Path


####### parameters ############

rasterPath="/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/Resistances/Hermine/resistances_glm16_experts_addition_clip.tif"

CScapeOutPath="/home/luvil/IE-OFEV/Hermine/CSCape_tiles_addition_model_glm16_experts"

tiles = pd.read_csv('/media/luvil/NAS_DEVELOPPEMENT/IE_OFEV/Pinch-points/Hermine/grille/grille_3000_hermine_t1.csv')

percentBuffer=100
tileSize=(3000,3000)
force_square=True
checkTotalSize=False
nProc=8

##############################


def writeINIfile(resTilePath,nodeTilePath, outPath, outName):
    return(
        f'[Options for advanced mode]\n'
f'ground_file_is_resistances = True\n'
f'source_file = (Browse for a current source file)\n'
f'remove_src_or_gnd = keepall\n'
f'ground_file = (Browse for a ground point file)\n'
f'use_unit_currents = True\n'
f'use_direct_grounds = False\n'
f'\n'
f'[Calculation options]\n'
f'low_memory_mode = False\n'
f'parallelize = True\n'
f'solver = cholmod\n'
f'print_timings = True\n'
f'max_parallel = 2\n'
f'\n'

f'[Options for pairwise and one-to-all and all-to-one modes]\n'
f'included_pairs_file = None\n'
f'use_included_pairs = False\n'
f'point_file ={nodeTilePath}\n'
f'\n'

f'[Output options]\n'
f'write_cum_cur_map_only = True\n'
f'log_transform_maps = False\n'
f'output_file = {os.path.join(outPath,outName)}\n'
f'write_max_cur_maps = False\n'
f'write_volt_maps = False\n'
f'set_null_currents_to_nodata = False\n'
f'set_null_voltages_to_nodata = False\n'
f'compress_grids = False\n'
f'write_cur_maps = True\n'
f'\n'

f'[Short circuit regions (aka polygons)]\n'
f'use_polygons = False\n'
f'polygon_file = None\n'
f'\n'

f'[Connection scheme for raster habitat data]\n'
f'connect_four_neighbors_only = True\n'
f'connect_using_avg_resistances = True\n'
f'\n'

f'[Habitat raster or graph]\n'
f'habitat_file = {resTilePath}\n'
f'habitat_map_is_resistances = True\n'
f'\n'

f'[Options for one-to-all and all-to-one modes]\n'
f'use_variable_source_strengths = False\n'
f'variable_source_file = None\n'
f'\n'
f'[Version]\n'
f'version = 4.0.5\n'
f'\n'
f'[Mask file]\n'
f'use_mask = False\n'
f'mask_file = None\n'
f'\n'
f'[Circuitscape mode]\n'
f'data_type = raster\n'
f'scenario = pairwise\n'
    )


def submitTile(rasterPath,rowStart, rowStop, colStart, colStop, padding, k, outputPath):
        with rio.open(rasterPath, 'r') as inputRasterStream:
                win=Window.from_slices((rowStart-padding, rowStop+padding), (colStart-padding, colStop+padding)) 
                win_transform = inputRasterStream.window_transform(win)
                rast=inputRasterStream.read(1, window=win)

                cellSize=win_transform[0]
                header={"ncols":rast.shape[1],
                        "nrows":rast.shape[0],
                        "xllcorner":int(win_transform[2]),
                        "yllcorner":int(win_transform[5]-(rast.shape[0]*cellSize)),
                        "cellSize":int(cellSize),
                        "NODATA_value":-9999}     

                northSouth=np.pad(rast[ 1:-1,:],((1,1), (0,0)), mode="constant", constant_values=1)
                northSouth[np.isnan(northSouth)]=-9999
                outFileResNS=os.path.join(outputPath, f"res_North_South_tile_{k}.txt")
                np.savetxt(outFileResNS,northSouth,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")


                eastWest=np.pad(rast[:,1:-1],((0,0), (1,1)), mode="constant", constant_values=1)
                eastWest[np.isnan(eastWest)]=-9999
                outFileResEW=os.path.join(outputPath,f"res_East_West_tile_{k}.txt")
                np.savetxt(outFileResEW,eastWest,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")

              
                nullMat=np.zeros(shape=rast.shape, dtype="int16")
                nullMat[nullMat==0]=-9999


                northSouth=np.pad(nullMat[ 1:-1,:],((1,1), (0,0)), mode="constant", constant_values=(1,2))
                eastWest=np.pad(nullMat[:,1:-1],((0,0), (1,1)), mode="constant", constant_values=(1,2))

                outFileNodeNS=os.path.join(outputPath, f"nodes_North_South_tile_{k}.txt")
                np.savetxt(outFileNodeNS,northSouth,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")
                outFileNodeEW=os.path.join(outputPath, f"nodes_East_West_tile_{k}.txt")
                np.savetxt(outFileNodeEW,eastWest,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")

                outFileCS_NS=f"t_{k}_ns.out"
                with open(os.path.join(outputPath, f"t{k}_ns.ini"), "x") as f:
                    f.write(writeINIfile(resTilePath=outFileResNS, nodeTilePath=outFileNodeNS,outPath=CScapeOutPath, outName=outFileCS_NS))
                outFileCS_EW=f"t_{k}_ew.out"
                with open(os.path.join(outputPath, f"t{k}_ew.ini"), "x") as f:
                    f.write(writeINIfile(resTilePath=outFileResEW, nodeTilePath=outFileNodeEW,outPath=CScapeOutPath, outName=outFileCS_EW))

                with open(os.path.join(outputPath, f"runTile{k}ns.jl"), "x") as f:
                        f.write('using Pkg\n')
                        f.write('Pkg.add("Circuitscape")\n')
                        f.write('Pkg.update()\n')
                        f.write('using Circuitscape\n')
                        f.write(f'compute("{os.path.join(outputPath, f"t{k}_ns.ini")}")')

                with open(os.path.join(outputPath, f"runTile{k}ew.jl"), "x") as f:
                        f.write('using Pkg\n')
                        f.write('Pkg.add("Circuitscape")\n')
                        f.write('Pkg.update()\n')
                        f.write('using Circuitscape\n')
                        f.write(f'compute("{os.path.join(outputPath, f"t{k}_ew.ini")}")')

        return(f" *** submitTile  for tile {k} done ****")
if __name__ == "__main__":    
        #pathName, fileName = os.path.split(outputPath)
        start=datetime.now()
        CScapeOutPath=Path(CScapeOutPath)
        tilesIndex=tiles["id"].to_list()

        if not os.path.exists(CScapeOutPath):
                print(""" ***** output directory for CircuitScape does not exists : please create it 
                (I could create it for you but I would not since output folder the size could be HUGE !!) ******
                """)
                sys.exit(0)
        with rio.open(rasterPath) as inp:
                expFactor=percentBuffer/100
                ncol=inp.width
                nrow=inp.height
                inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
                cellSize=inp_transform[0]
                
                print("#"*10)
                print("input raster features : ")
                print(f"{nrow} row X {ncol} columns ")
                print(f"cell size : {cellSize}")
                print(f"'nodata' code : {inp.nodata}")

                ### input parameters
                padding=int(tileSize[0]*expFactor) ### the number of cells to add on the contours
                print(f"buffer size : {padding} cells")
                totalSize=np.power(padding*2+tileSize[0],2)
                print(f"total size of a tile is {padding*2+tileSize[0]} rows X {padding*2+tileSize[1]} columns")
                print(f"total number of pixels for tile + buffer size : {totalSize} cells")
                if checkTotalSize:
                        assert  totalSize < 10e8, "the tile Size + buffer is greater than 100'000'000 cells : reduce the tile size or the buffer percentage"
                
                if((ncol%tileSize[0]!=0) or (nrow%tileSize[1]!=0)):
                       print("the tileSize does not fully divide the raster")
                       if force_square:
                                print("you can choose among the followest : ")
                                print(np.gcd(ncol, nrow))
                                sys.exit(0)  
                              
                              
                                   
                              

                nbSplitX=ncol//(tileSize[1])
                nbSplitY=nrow//(tileSize[0])
                print("number of tiles :")
                print(f"{nbSplitX} tiles per row X {nbSplitY} tiles per column")
                
                #assert ncol%tileSize[0]==0,"the remainder of the divison for the X axis is not zero : the "
                splitCoordsXStart=[ i[0] for i in np.array_split(np.arange(ncol),nbSplitX)]
                splitCoordsXStop=[ i[-1] for i in np.array_split(np.arange(ncol),nbSplitX)]
                splitCoordsYStart=[ i[0] for i in np.array_split(np.arange(nrow),nbSplitY)]
                splitCoordsYStop=[ i[-1] for i in np.array_split(np.arange(nrow),nbSplitY)]
                k=0

        with ProcessPoolExecutor(nProc) as executor:
                futures=[]
                for i in range(len(splitCoordsXStart)):
                        for j in range(len(splitCoordsYStart)):
                                k=k+1
                                colStart=splitCoordsXStart[i]
                                colStop=splitCoordsXStop[i]+1
                                rowStart=splitCoordsYStart[j]
                                rowStop=splitCoordsYStop[j]+1

                                if (rowStart-padding<0 or colStart-padding<0 or colStop+padding>ncol or rowStop+padding>nrow):
                                    print(f"the tile {k} is a border tile : no computation required")
                                    continue
                                
                               
                        #   elif rowStart-padding<0: ## these are for the tiles that are on top (apart from the upper-left tile)
                        #          win=Window.from_slices((rowStart, rowStop+padding), (colStart-padding, colStop+padding))
                        #  elif colStart-padding<0: ## these are for the tiles that are on the left (apart from the upper-left tile)
                        #         win=Window.from_slices((rowStart-padding, rowStop+padding), (colStart, colStop+padding))
                                else:
                                        
                                        if k not in tilesIndex:
                                            print(f' tile {k} is not in the tile indexes you provided :: skipping it')
                                            continue
                                        else:
                                            futures.append(executor.submit(submitTile, rasterPath, rowStart, rowStop,colStart, colStop, padding, k , CScapeOutPath))
                for r in as_completed(futures):
                        print(r.result())
        print(" >>>> SCRIPT COMPLETED")
        print(f" total duration : {datetime.now()-start}")
        print(f" k = {k}")
