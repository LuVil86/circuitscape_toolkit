    
import rasterio as rio
import numpy as np
import scipy.ndimage as ndi
inputRaster= "/media/luvil/T7_ROUGE/tuile_1/ZN_tile_155.txt"
outSource = "/media/luvil/T7_ROUGE/tuile_1/source_patch_edges_tile_155.txt"
#outSource = "/media/luvil/T7_ROUGE/tuile_1/only_patch_edges.tif"


with rio.open(inputRaster) as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
        cellSize=inp_transform[0]
    
        matrix=inp.read(1)
    

        
        maskBin = np.where(matrix==-9999, 0,1)
        ### le masque de l'érosion donne les endroits où l'algorithme doit opérer (il évitera les autres)
        ## -> on demande de faire l'érosion que sur les chemins
        
        #finMask = ndi.binary_erosion(maskBoth,iterations=5,mask=maskRoad, brute_force=False)
        erosion = ndi.binary_erosion(maskBin, iterations=1)
        add = np.add(erosion, maskBin)
        diff = np.array([[elem in [0,2] for elem in row] for row in add], dtype=np.uint8) 
        #finArray = np.where(diff, -9999, 1)
        finArray=np.where(diff,-9999, matrix)
        #finArray = np.where(erosion ==2, 1, -9999)

        header={"ncols":finArray.shape[1],
        "nrows":finArray.shape[0],
        "xllcorner":int(inp_transform[2]),
        "yllcorner":int(inp_transform[5]-(finArray.shape[0]*cellSize)),
        "cellSize":int(cellSize),
        "NODATA_value":-9999}   

        np.savetxt(outSource,finArray,delimiter='\t', header="\n".join([f"{i}\t{header[i]}" for i in header]),comments="", encoding="utf-8", fmt="%d")



