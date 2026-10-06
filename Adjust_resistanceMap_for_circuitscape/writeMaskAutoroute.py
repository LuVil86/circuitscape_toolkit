    
import rasterio as rio
import numpy as np
import scipy.ndimage as ndi
inputRaster= "/media/luvil/T7_ROUGE/tuile_1/habitatMap_tile_155.tif"
inputRes = "/media/luvil/T7_ROUGE/tuile_1/res_East_West_tile_155.txt"
outGround = "/media/luvil/T7_ROUGE/tuile_1/ground_autoroute_tile_155.txt"
outMask = "/media/luvil/T7_ROUGE/tuile_1/mask_autoroute_tile_155.txt"
outSource = "/media/luvil/T7_ROUGE/tuile_1/source_autoroute_tile_155.txt"
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
        
        maskBin = np.where(matrix==23, 1,0)
        maskAuto =np.where(matrix==23, -9999,1)
        maskSource =np.where(matrix==23, -1000,-9999)
        maskImp = np.where(maskBin == 1, -9999, resMap)
        ### le masque de l'érosion donne les endroits où l'algorithme doit opérer (il évitera les autres)
        ## -> on demande de faire l'érosion que sur les chemins
        
        #finMask = ndi.binary_erosion(maskBoth,iterations=5,mask=maskRoad, brute_force=False)
        erosion = ndi.binary_erosion(maskBin, iterations=1)
        add = np.add(erosion, maskBin)
        diff = np.array([[elem in [0,2] for elem in row] for row in add], dtype=np.uint8) 
        finArray = np.where(diff, -9999, 1)

        header={"ncols":finArray.shape[1],
        "nrows":finArray.shape[0],
        "xllcorner":int(inp_transform[2]),
        "yllcorner":int(inp_transform[5]-(finArray.shape[0]*cellSize)),
        "cellSize":int(cellSize),
        "NODATA_value":-9999}   
        np.savetxt(outGround,finArray,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")
        np.savetxt(outMask,maskAuto,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")
        np.savetxt(outSource,maskSource,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")
        np.savetxt(outImpassable,maskImp,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")

