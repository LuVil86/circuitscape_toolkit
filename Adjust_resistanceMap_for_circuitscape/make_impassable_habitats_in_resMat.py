    
import rasterio as rio
import numpy as np
import scipy.ndimage as ndi
inputRaster= "/media/luvil/T7_ROUGE/tuile_1/habitatMap_tile_155.tif"
inputRes = "/media/luvil/T7_ROUGE/tuile_1/res_East_West_tile_155.txt"
outImpassable = "/media/luvil/T7_ROUGE/tuile_1/resMat_autoroute_impassable_tile_155.txt"
with rio.open(inputRaster) as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
        cellSize=inp_transform[0]
    
        matrix=inp.read(1)
    
with rio.open(inputRes) as res:

        resMap = res.read(1)
        
        maskBin = np.array([[elem in [23] for elem in row] for row in matrix], dtype=np.uint8) 
        finArray = np.where(maskBin, -9999, resMap)
        ### le masque de l'érosion donne les endroits où l'algorithme doit opérer (il évitera les autres)
        ## -> on demande de faire l'érosion que sur les chemins
        
        header={"ncols":finArray.shape[1],
        "nrows":finArray.shape[0],
        "xllcorner":int(inp_transform[2]),
        "yllcorner":int(inp_transform[5]-(finArray.shape[0]*cellSize)),
        "cellSize":int(cellSize),
        "NODATA_value":-9999}   
        np.savetxt(outImpassable,finArray,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")

