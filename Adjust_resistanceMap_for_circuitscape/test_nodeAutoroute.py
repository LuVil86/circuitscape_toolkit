
import rasterio as rio
import numpy as np
import scipy.ndimage as ndi
inputRaster= "/media/luvil/T7_ROUGE/tuile_1/habitatMap_tile_155.tif"
inputNodes="/media/luvil/T7_ROUGE/tuile_1/nodes_East_West_tile_155.txt"
outNodes = "/media/luvil/T7_ROUGE/tuile_1/nodes_autoroute_tile_155.txt"

with rio.open(inputRaster) as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform 
        cellSize=inp_transform[0]
    
        matrix=inp.read(1)
        maskBin = np.array([[elem == 23 for elem in row] for row in matrix], dtype=np.uint8) 



with rio.open(inputNodes) as inp:
        nodes = inp.read(1)
        finArray = np.where(maskBin==True, 3, nodes)
        header={"ncols":finArray.shape[1],
        "nrows":finArray.shape[0],
        "xllcorner":int(inp_transform[2]),
        "yllcorner":int(inp_transform[5]-(finArray.shape[0]*cellSize)),
        "cellSize":int(cellSize),
        "NODATA_value":-9999}   
        np.savetxt(outNodes,finArray,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")
      