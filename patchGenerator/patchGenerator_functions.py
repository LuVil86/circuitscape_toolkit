import numpy as np
import pandas as pd
from osgeo import gdal, osr
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from patchify import patchify



### parse matrix tiff ###
def parse_matrix(path, file, nb_mat):    #nb_mat is the number of sub matrix
    ds = gdal.Open(path + file, gdal.GA_ReadOnly)
    dictcoor = {}
    nrow = ds.RasterYSize
    ncol = ds.RasterXSize
    print("Original matrix size: ", (nrow, ncol))
    d = int(np.sqrt(nb_mat))
    stepx = int(ncol/d)
    stepy = int(nrow/d)
    xmin = 0
    ymin = 0
    xsteps = [xmin + stepx * i for i in range(d)]
    ysteps = [ymin + stepy * i for i in range(d)]
    for i in range(d):
        for j in range(d):
            dictcoor[(i,j)] = [xsteps[j], ysteps[i], stepx, stepy]
    return dictcoor


### from raster to matrix ###
def raster_to_matrix(path,file,x,y,stepx,stepy):                                  #previously "dame_matriz"
    ds = gdal.Open(path + file, gdal.GA_ReadOnly)
    matrice = np.array(ds.GetRasterBand(1).ReadAsArray(x, y, stepx,stepy), dtype=np.int16)
    return matrice

### from raster to matrix  OLD oNE###
def raster_to_matrix_old(path,file):
    ds = gdal.Open(path+file, gdal.GA_ReadOnly)
    matrice = np.array(ds.GetRasterBand(1).ReadAsArray())
    return matrice


### extract the header of the txt file ###
def header_txt(path,file,separ):                     #separ = '\t' or "," or " "        previously "dame_cabecera"
    f = open(path + file,'r')                        #open txt file, 'r' = read
    dictio = {}                                         #create empty dictionary
    for line in f:                                      #read each line of the text file
        aux = line.split(separ)                         #separe (according to the chosen seaprator) the line in a list
        if len(aux) < 20 :                              #if the list is shorter than 20 (why not just 3?)
            dictio.update({aux[0]:int(aux[-1][:-1])})   #[-1][:-1] : takes the last element of the list. In this element take everything except the last character (because it's a line break)
        else :
            break
    f.close()
    return dictio                                       #return example : {'ncols': 6940, 'nrows': 10940, 'xllcorner': 2655300, 'yllcorner': 1217800, 'cellsize': 5, 'NODATA_value': -9999}

#load matrix (from tkt file ?)
def matrix_nodes(path,file,separ):                      #previously "matriz_nodos"
    with open(path+file, 'r') as f:
        l = [[int(num) for num in line.split(separ)] for line in f if line[0].isdigit() ]   #list Comprehension. load the file into a matrix (array) if the first element of the line is a digit (number)
    return np.array(l)

def show_image(matrix):                                 #previously "muestra_image"
    plt.imshow(matrix)
    plt.show()

#previously "elimino_caminos"                     ### remove the pathways when they are next to (to define) a milieu of interest (replace the code "pathway" by the code "milieu") ###
def clean_pathway(matrix,pathway, milieu, filter_size = 3) :                            #matrix (matrix) = matrix of the habitats | pathway (integer) = code for pathway in the matrix | milieu (integer) = code for the milieu (habitat) in the matrix
    fs = filter_size
    add_mat = int((fs-1) / 2)                                                           #add_mat = 1 (size of the edge added arond the matrix)
    matrix_expansa = np.zeros((matrix.shape[0]+add_mat*2,matrix.shape[1]+add_mat*2))    #create empty matrix (full of 0). nb line = nb line of "matrix"+ 2. nb column = nb columns of "matrix"+ 2 (need *2 to have a  line more at left/right and up/down)
    matrix_expansa[add_mat:-add_mat,add_mat:-add_mat] = matrix                          #inlude "matrix" inside matrix_expensa. The size of the buffer around the matrix depend of "add_mat"   Example  [[0. 0. 0. 0. 0.]            
    mat_pathway = matrix_expansa == pathway                                             #extract the pathways from matrix_expensa (where matrix_expensa = pathway code => True)                          [0. 1. 1. 1. 0.]
    mat_milieu = matrix_expansa == milieu                                                #extract the milieu from matrix_expensa (where matrix_expensa = milieu code => True)                            [0. 1. 1. 1. 0.] 
    for i in range(add_mat,matrix_expansa.shape[0]-add_mat) :                           #shape[0] = row                                                                                                  [0. 1. 1. 1. 0.]
        for j in range(add_mat,matrix_expansa.shape[1]-add_mat):                        #shape[1] = column                                                                                               [0. 0. 0. 0. 0.]]
            aux = matrix_expansa[i-add_mat:i+add_mat+1, j-add_mat:j+add_mat+1]          #define the area around the point (+1 to have a 3x3 matrix with element i,j in the center)                                                                         
            if mat_pathway[i,j] and milieu in aux:                                      #if pathway next to milieu
                mat_milieu[i,j] = True                                                       #in the milieu matrix, what was count as pathway (False) is now count as milieu (True)
    matrix_expansa[mat_milieu] = milieu                                                  #in the matrix_expensa, where mat_milieu = True, replace by code milieu

    return matrix_expansa[add_mat:-add_mat,add_mat:-add_mat]                            #remove the edge


#previously "elimino_caminos_multi"                 ### remove the pathways when they are next to (to define) a milieu of interest (replace the code "pathway" by the code "milieu") ###
def clean_pathway_multi(data) :                                                         #data  = list
    matrix = data[0]                                                                    #matrix (matrix) = matrix of the habitats
    pathway = data[1]                                                                   #pathway (integer) = code for pathway in the matrix
    milieu = data[2]                                                                    #milieu (integer) = code for the milieu (habitat) in the matrix
    filter_size = data[3]
    fs = filter_size
    add_mat = int((fs-1) / 2)                                                           #add_mat = 1 (size of the edge added arond the matrix)
    matrix_expansa = np.zeros((matrix.shape[0]+add_mat*2,matrix.shape[1]+add_mat*2))    #create empty matrix (full of 0). nb line = nb line of "matrix"+ 2. nb column = nb columns of "matrix"+ 2 (need *2 to have a  line more at left/right and up/down)
    matrix_expansa[add_mat:-add_mat,add_mat:-add_mat] = matrix                          #inlude "matrix" inside matrix_expensa. The size of the buffer around the matrix depend of "add_mat"   Example  [[0. 0. 0. 0. 0.]            
    mat_pathway = matrix_expansa == pathway                                             #extract the pathways from matrix_expensa (where matrix_expensa = pathway code => True)                          [0. 1. 1. 1. 0.]
    mat_milieu = matrix_expansa == milieu                                                #extract the milieu from matrix_expensa (where matrix_expensa = milieu code => True)                            [0. 1. 1. 1. 0.] 
    for i in range(add_mat,matrix_expansa.shape[0]-add_mat) :                           #shape[0] = row. for in range (1,nrow-1)                                                                         [0. 1. 1. 1. 0.]
        for j in range(add_mat,matrix_expansa.shape[1]-add_mat):                        #shape[1] = column for j in range (1,ncol-1)                                                                     [0. 0. 0. 0. 0.]]
            aux = matrix_expansa[i-add_mat:i+add_mat+1, j-add_mat:j+add_mat+1]          #define the area around the point (+1 to have a 3x3 matrix with element i,j in the center)                                                                         
            if mat_pathway[i,j] and milieu in aux:                                      #if pathway next to milieu
                mat_milieu[i,j] = True                                                       #in the milieu matrix, what was count as pathway (False) is now count as milieu (True)
    matrix_expansa[mat_milieu] = milieu                                                  #in the matrix_expensa, where mat_milieu = True, replace by code milieu

    return matrix_expansa[add_mat:-add_mat,add_mat:-add_mat]                            #remove the edge

#previously  "graba_ascii"                          ### save matrix as ascii (txt with header) ###
def write_ascii(name,header,matrix,separ):
    g = open(name,'w')
    for k,v in header.items():
        g.write(separ.join([k,str(v)]))
        g.write('\n')
    matrix = np.array(matrix,dtype=int)
    for i in range(matrix.shape[0]):
        g.write(separ.join(np.array(matrix[i,:],dtype=str)))
        g.write('\n')
    g.close()


def write_geotiff(new_array, output_file, original_file):
    ds = gdal.Open(original_file, gdal.GA_ReadOnly)
    band = ds.GetRasterBand(1)
    geotransform = ds.GetGeoTransform()
    wkt = ds.GetProjection()

    # Create gtif file
    driver = gdal.GetDriverByName("GTiff")

    dst_ds = driver.Create(output_file,
                        band.XSize,
                        band.YSize,
                        1,gdal.GDT_Int32)   #gdal.GDT_Int16

    #writting output raster
    dst_ds.GetRasterBand(1).WriteArray(new_array)
    #setting nodata value
    dst_ds.GetRasterBand(1).SetNoDataValue(-999)
    #setting extension of output raster
    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    dst_ds.SetGeoTransform(geotransform)
    # setting spatial reference of output raster
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    dst_ds.SetProjection( srs.ExportToWkt() )
    #Close output raster dataset

    ds = None
    dst_ds = None




### use in "clusteriz" to find the number of the cluster ###
def find_clust_num_f(dove):     #dove = list len(dove)=2 (the element upper and the one of the left)
    dove.sort()                 #sort in ascending order
    if dove.count(0) == 2:      #if the list countains 2 zeros
        return [0]                  #array with one element : "0"
    elif dove.count(0) == 1 :   #if the list countains 1 zero
        return [dove[1]]            #takes the second element of the list (bcse the first one is 0)
    else:                       # else = the list doesn't countain any 0
        return dove                 #return all the list dove
    



#previously "clusterizo"                                    ### create patches when pixel of the same milieu are adjacent (up and left but not in diagonal) ###
def clusteriz(matrix) :  #, filter_size = 3                                      #input matrix = matrix milieu (only one milieu)
    clust_num = 0                                                               #number of the cluster
    fs =  3            #filter_size  =                                              
    add_mat = int((fs-1) / 2)                                                   #add_mat = 1 (size of the edge added arond the matrix)
    m_big = np.zeros((matrix.shape[0]+add_mat*2,matrix.shape[1]+add_mat*2),dtype=int)     #create empty matrix (full of 0). nb line = nb line of "matrix"+ 2. nb column = nb columns of "matrix"+ 2 (need *2 to have a  line more at left/right and up/down)
    m_big[add_mat:-add_mat,add_mat:-add_mat] = matrix                           #inlude "matrix" inside m_big. The size of the buffer around the matrix depend of "add_mat"   
    m_cluster = np.zeros(m_big.shape, dtype=int)                                           #create array full of 0 with the shape of m_big (m_cluster = matrix of the cluster (will be filled in))
    for i in np.arange(1,m_big.shape[0]-1):                                     #go through the matrix m_big    #shape[0] : row             #why 1, -1 and not add_mat, -add_mat ? (like for clean pathway)
        for j in np.arange(1, m_big.shape[1]-1):                                                                #shape[1] : column
            if m_big[i, j] != 0:                                                #if the element (at i, j) in matrix m_big is not equal to 0
                il = [m_cluster[i - 1, j], m_cluster[i, j - 1]]                        # il = list with two elements (the element upper and the one of the left (in the cluster matrix) of m_big[i,j])
                cuales = find_clust_num_f(il)                                               # if il = [0,0] => find_clust_num_f(il) = [0] ; if il = [1,0] => find_clust_num_f(il) = [1] ; if il = [1,1] => find_clust_num_f(il) = [1,1]
                if len(cuales) == 2:
                    if cuales[0] != cuales[1] :                                 #for example cuales = [1,2] (two clusters around)
                        # print(cuales,cuales[0],cuales[1])                        
                        m_cluster[i, j] = cuales[0]                                 #in matrix cluster (at same position than in matrix big), the element takes the number of the cluster (here n°1)
                        m_cluster[m_cluster == cuales[1]] = cuales[0]               #in matrix cluster, where the element equal to "2" replace by "1" because if adjacent to n°1 means they are part of the same cluster. example where it's necessary      0 0 0 0 0 0    
                    else:                                                                                                           #                                                                                                                       0 0 0 1 1 0
                        m_cluster[i, j] = cuales[0]                             #same cluster up and on the left                    #                                                                                                                       0 2 2 x   0   x will become 1 so the two 2 needs also to become 1
                else:
                    if cuales[0] == 0:                                          #if no cluster around (at the begining, cluster matrix full of 0, so il = [0,0] => cuales = [0])
                        clust_num += 1                                                      #count number of cluster, create a new one                                              
                        m_cluster[i, j] = clust_num                             
                    else:                                                       #if just one cluster around (cuales[0] = n° of the cluster)
                        m_cluster[i, j] = cuales[0]                         
    print(len(np.unique(m_cluster)))                                      #len(np.unique(m_cluster)) gives the number of cluster
    return m_cluster[add_mat:-add_mat,add_mat:-add_mat]                                                            #return matrix cluster


#previously "clusterizo_multi"                               ### create patches when pixel of the same milieu are adjacent (up and left but not in diagonal) ###
def clusteriz_multi(data) :
    start = datetime.now()
    matrix =  data[0]                                                           #input matrix = matrix milieu
    filter_size = data[1]
    milieu = data[2]                                                            #milieu of interest (to be clusterized)
    #clust_num = 1                                                               #number of the cluster
    fs = filter_size
    add_mat = int((fs-1) / 2)                                                   #add_mat = 1 (size of the edge added arond the matrix)
    m_big = np.zeros((matrix.shape[0]+add_mat*2,matrix.shape[1]+add_mat*2), dtype=int)     #create empty matrix (full of 0). nb line = nb line of "matrix"+ 2. nb column = nb columns of "matrix"+ 2 (need *2 to have a  line more at left/right and up/down)
    m_big[add_mat:-add_mat,add_mat:-add_mat] = matrix                           #inlude "matrix" inside m_big. The size of the buffer around the matrix depend of "add_mat"   
    m_cluster = np.zeros(m_big.shape, dtype=int)                                           #create array full of 0 with the shape of m_big (m_cluster = matrix of the cluster (will be filled in))
    for i in np.arange(1,m_big.shape[0]-1):                                     #go through the matrix m_big    #shape[0] : row             #why 1, -1 and not add_mat, -add_mat ? (like for clean pathway)
        for j in np.arange(1, m_big.shape[1]-1):                                                                #shape[1] : column
            if m_big[i, j] != 0:                                                #if the element (at i, j) in matrix m_big is not equal to 0
                il = [m_cluster[i - 1, j], m_cluster[i, j - 1]]                     # il = list with two elements (the element upper and the one of the left (in the cluster matrix) of m_big[i,j])
                cuales = find_clust_num_f(il)                                       # if il = [0,0] => find_clust_num_f(il) = [0] ; if il = [1,0] => find_clust_num_f(il) = [1],  if il = [1,1] => find_clust_num_f(il) = [1,1]
                if len(cuales) == 2:
                    if cuales[0] != cuales[1] :                                 #for example cuales = [1,2] (two clusters around)
                        # print(cuales,cuales[0],cuales[1])                        
                        m_cluster[i, j] = cuales[0]                                 #in matrix cluster (at same position than in matrix big), the element takes the number of the cluster (here n°1)
                        m_cluster[m_cluster == cuales[1]] = cuales[0]               #in matrix cluster, where the element equal to "2" replace by "1" because if adjacent to n°1 means they are part of the same cluster. example where it's necessary      0 0 0 0 0 0    
                    else:                                                                                                           #                                                                                                                       0 0 0 1 1 0
                        m_cluster[i, j] = cuales[0]                             #same cluster up and on the left                    #                                                                                                                       0 2 2 x   0   x will become 1 so the two 2 needs also to become 1
                else:
                    if cuales[0] == 0:                                          #if no cluster around (at the begining, cluster matrix full of 0, so il = [0,0] => cuales = [0])
                        clust_num = m_cluster.max()+1                                         #there was a pb with the clust_num, I tried to take the max value of the new cluster matrix to find the next cluster number
                        m_cluster[i, j] = clust_num                              
                    else:                                                       #if just one cluster around (cuales[0] = n° of the cluster)
                        m_cluster[i, j] = cuales[0]                                 #takes the value of the cluster around
                                                             
    print('The clustering of the milieu {} is done. It has {} patches and it tooks {}'.format(milieu,len(np.unique(m_cluster))-1,datetime.now()-start))
    return m_cluster[add_mat:-add_mat,add_mat:-add_mat]

#TEST NOT COMPLETED WITH SLIDING WINDOWS                              ### create patches when pixel of the same milieu are adjacent (up and left but not in diagonal) ###
def clusteriz_multi_wind(data) :
    start = datetime.now()
    matrix =  data[0]                                                           #input matrix = matrix milieu
    filter_size = data[1]
    milieu = data[2]                                                            #milieu of interest (to be clusterized)
    wind_size = data[3]
    step = data[4]
    fs = filter_size
    add_mat = int((fs-1) / 2)                                                   #add_mat = 1 (size of the edge added arond the matrix)
    m_big = np.zeros((matrix.shape[0]+add_mat*2,matrix.shape[1]+add_mat*2), dtype=int)     #create empty matrix (full of 0). nb line = nb line of "matrix"+ 2. nb column = nb columns of "matrix"+ 2 (need *2 to have a  line more at left/right and up/down)
    m_big[add_mat:-add_mat,add_mat:-add_mat] = matrix                           #inlude "matrix" inside m_big. The size of the buffer around the matrix depend of "add_mat"   
    m_cluster = np.zeros(m_big.shape, dtype=int)                                           #create array full of 0 with the shape of m_big (m_cluster = matrix of the cluster (will be filled in))
    patches=patchify(m_big, (wind_size[0],wind_size[1]), step=step)             #sliding windows
    n_mat_row = patches.shape[0] #number of matrix in row
    n_mat_col = patches.shape[1] #number of matrix in col

    for x in range(n_mat_row):
        for y in range(n_mat_col):
            sub_mat = patches[x,y]
            gapi = step*x
            gapj = step*y
    
            for i in range(1,sub_mat.shape[0]-1):                                     #go through the matrix m_big    #shape[0] : row             #why 1, -1 and not add_mat, -add_mat ? (like for clean pathway)
                for j in np.arange(1, sub_mat.shape[1]-1):                                                                #shape[1] : column
                    if sub_mat[i, j] != 0:                                                #if the element (at i, j) in matrix m_big is not equal to 0
                        il = [m_cluster[i+gapi - 1, j+gapj], m_cluster[i+gapi, j+gapj - 1],m_cluster[i+gapi - 1, j+gapj-1]]                     # il = list with two elements (the element upper and the one of the left (in the cluster matrix) of m_big[i,j])
                        cuales = find_clust_num_diag(il)                                       # if il = [0,0] => find_clust_num_f(il) = [0] ; if il = [1,0] => find_clust_num_f(il) = [1],  if il = [1,1] => find_clust_num_f(il) = [1,1]
                        if len(cuales) == 3:
                            if cuales[0] != cuales[1] :                                 #for example cuales = [1,2] (two clusters around)
                                # print(cuales,cuales[0],cuales[1])                        
                                m_cluster[i+gapi, j+gapj] = cuales[0]                                 #in matrix cluster (at same position than in matrix big), the element takes the number of the cluster (here n°1)
                                m_cluster[m_cluster == cuales[1]] = cuales[0]               #in matrix cluster, where the element equal to "2" replace by "1" because if adjacent to n°1 means they are part of the same cluster. example where it's necessary      0 0 0 0 0 0    
                            else:                                                                                                           #                                                                                                                       0 0 0 1 1 0
                                m_cluster[i+gapi, j+gapj] = cuales[0]                             #same cluster up and on the left                    #                                                                                                                       0 2 2 x   0   x will become 1 so the two 2 needs also to become 1
                        else:
                            if cuales[0] == 0:                                          #if no cluster around (at the begining, cluster matrix full of 0, so il = [0,0] => cuales = [0])
                                clust_num = m_cluster.max()+1                                         #there was a pb with the clust_num, I tried to take the max value of the new cluster matrix to find the next cluster number
                                m_cluster[i+gapi, j+gapj] = clust_num                              
                            else:                                                       #if just one cluster around (cuales[0] = n° of the cluster)
                                m_cluster[i+gapi, j+gapj] = cuales[0]                                 #takes the value of the cluster around
                                                                
    print('The clustering of the milieu {} is done. It has {} patches and it tooks {}'.format(milieu,len(np.unique(m_cluster))-1,datetime.now()-start))
    return m_cluster[add_mat:-add_mat,add_mat:-add_mat]



#previously "lee_areas_multi"                                     ### check if the patch is big enough to be taken into consideration ###
def lee_areas_multi(data):
    pixel_size = 5
    file = data[0]                                                         #csv file with the patches for one milieu (ex: file = 'medio_2.csv' (need absolutely to be in this format bcse "milieu = int(file[6:-4])"")
    path = data[1]                                                         #working directory
    areas_lim = data[2]                                                    #dictionary with code milieu / minimum patch size (ex {'2':30 ,'4': 30,'15':30})

    print('=' * 60)
    milieu = int(file[6:-4])                                               #takes the number between "_" and ".csv" in the file name. Here it's "2"
    print('milieu: ', milieu)                                                
    print(file)

    matrice = np.loadtxt(path + file, delimiter=',')                       #matrice of the patches for one milieu. it's "loadtxt" but file = ".csv" with "," as delimiter
    patches = np.unique(matrice)[np.unique(matrice) != 0]                   #all the value for patch

    print('=' * 60)                                                                   #if {:02d} : fill with one zero before if the number is one len 1. ex : '{:02d}'.format(3) = 03
    print('There are {1:2d} patches of milieu {0:2d}'.format(milieu, len(patches)))     #{1:2d} : element 1 (len(patches))  2 : with minimum size two (can be 'blank'number), formated in base 10 (d)
    print('-' * 60)
    # print('The total extent of milieu {:2d} is {:3d}'.format(milieu, np.count_nonzero(m_big == milieu)))
    area = pd.DataFrame(columns=('patch', 'area'))                                    #create new df with 2 columns called 'patch' and 'area'
    count = 0
    dictio = {}
    for i, el in enumerate(patches):                                                     #enumerate() simplify Looping With Counters. i = 0 , el = n° of the patch
        dove = matrice == el                                                            #dove = matrice with only one patch ("true" in place of the patch and "false" everywhere else)
        super = np.count_nonzero(dove) * pixel_size ** 2 * 1e-4                         #np.count_nonzero(dove): count number of "true"
        if super > areas_lim[str(milieu)]:                                              #if the size of the patch is biger that the limit size for this milieu
            print('patch n°:{} from milieu {} is {}° that fulfills the size condition, number of the patch:{}, size of the patch: {:3.3f} '.format(i,milieu,count,el,super))
            dictio[count] = {'patch': el, 'area': super}                                    #ex: {0: {'patch': 2, 'area': 234}, 1: {'patch': 5, 'area': 456}} 
            count += 1

    if not dictio:                                                                     #if dictio is empty
        print('There is no patch big enough for milieu:', milieu)
    else:                                                                              #if dictio not empty
        area = pd.DataFrame.from_dict(dictio.values(), orient='columns')                    #fill the df "area" with the value in dictio
        area = area.sort_values(by='area', ascending=False)                                 #sort by size (descending)
        print('=' * 60)
        print('How many patches are big enough', area.shape[0])                                               #nb of rows = nb of patches
        print('=' * 60)
        area.to_csv(path + file[6:-4] + '_area_info(newtry).csv', index=False)                 #file[6:-4] gives the number of the milieu
        print('End milieu ', milieu)







