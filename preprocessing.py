# Import libraries
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from osgeo import gdal, ogr, osr
gdal.UseExceptions() # So the errors are exceptions, and not just print errors

################# FUNCTIONS #################
# Rasterize counties and match to BUA 1975 file
def rasterize_counties(input_vfile, output_file, template_file):
    gdf=gpd.read_file(input_vfile)
    gdf['GEOID_INT']=gdf['GEOID'].map(int)
    gdf.to_file(input_vfile)
    
    ras_ds_template = gdal.Open(template_file) # use BUA1975_file as the template file for counties
    vec_ds = ogr.Open(input_vfile)
    lyr = vec_ds.GetLayer()    
    geot = ras_ds_template.GetGeoTransform()
    
    chn_ras_ds = gdal.GetDriverByName("GTiff").Create(output_file, 
                                                      ras_ds_template.RasterXSize, 
                                                      ras_ds_template.RasterYSize, 
                                                      1, 
                                                      gdal.GDT_UInt16, 
                                                      options = ['TFW=YES', 'COMPRESS=LZW'])
    chn_ras_ds.SetGeoTransform(geot)
    chn_ras_ds.FlushCache()
    
    chn_ras_ds.GetRasterBand(1).SetNoDataValue(0.0) 
    chn_ras_ds.SetSpatialRef(ras_ds_template.GetSpatialRef())
    gdal.RasterizeLayer(chn_ras_ds, [1], lyr, options=['ATTRIBUTE=GEOID_INT'])
    chn_ras_ds = None

    del chn_ras_ds
    del vec_ds
    
    testarr = gdal.Open(output_file).ReadAsArray()
    print(np.unique(testarr))
    del testarr

# Reproject counties to Mollweide projection, to be used for GHS preprocessing
def reproject_vector(input_file, in_epsg, output_file, out_epsg):
    driver = ogr.GetDriverByName('ESRI Shapefile')

    # input SpatialReference
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.SetFromUserInput(in_epsg)

    # output SpatialReference
    outSpatialRef = osr.SpatialReference()
    outSpatialRef.SetFromUserInput(out_epsg)

    # create the CoordinateTransformation
    coordTrans = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

    # get the input layer
    inDataSet = driver.Open(input_file)
    inLayer = inDataSet.GetLayer()

    # create the output layer
    if os.path.exists(output_file):
        driver.DeleteDataSource(output_file)
    outDataSet = driver.CreateDataSource(output_file)
    outLayer = outDataSet.CreateLayer("reproject", outSpatialRef, geom_type=ogr.wkbMultiPolygon)

    # add fields
    inLayerDefn = inLayer.GetLayerDefn()
    for i in range(0, inLayerDefn.GetFieldCount()):
        fieldDefn = inLayerDefn.GetFieldDefn(i)
        outLayer.CreateField(fieldDefn)

    # get the output layer's feature definition
    outLayerDefn = outLayer.GetLayerDefn()

    # loop through the input features
    inFeature = inLayer.GetNextFeature()
    while inFeature:
        # get the input geometry
        geom = inFeature.GetGeometryRef()
        # reproject the geometry
        geom.Transform(coordTrans)
        # create a new feature
        outFeature = ogr.Feature(outLayerDefn)
        # set the geometry and attribute
        outFeature.SetGeometry(geom)
        for i in range(0, outLayerDefn.GetFieldCount()):
            outFeature.SetField(outLayerDefn.GetFieldDefn(i).GetNameRef(), inFeature.GetField(i))
        # add the feature to the shapefile
        outLayer.CreateFeature(outFeature)
        # dereference the features and get the next input feature
        outFeature = None
        inFeature = inLayer.GetNextFeature()

    # Save and close the shapefiles
    inDataSet = None
    outDataSet = None

# Extract GHS data with buffer
def extract_GHS(input_file, output_file, cutline_file):
    # Assign boundary cutline
    warpOptions = gdal.WarpOptions(format = 'GTiff',
                                   cutlineDSName = cutline_file,
                                   cropToCutline=True,
                                   outputType = gdal.GDT_UInt16,
                                   dstNodata = 65535,
                                   creationOptions = ['TFW=YES', 'COMPRESS=LZW'],
                                   resampleAlg = 'near')
    # Extract and save
    gdal.Warp(output_file, input_file, options=warpOptions)
    

# Function to oversample GHS layer in 100 m resolution to 10 m resolution
def oversample_GHS(input_file, output_file):
    # Open raster and get band
    in_ds = gdal.Open(input_file)
    in_band = in_ds.GetRasterBand(1)
    
    # Multiply output size by 10
    out_rows = in_band.YSize * 10
    out_columns = in_band.XSize * 10

    # Create new data source (raster)
    gtiff_driver = gdal.GetDriverByName('GTiff')
    out_ds = gtiff_driver.Create(output_file, out_columns, out_rows, 1, gdal.GDT_Float32, options = ['TFW=YES', 'COMPRESS=LZW'] ) 
    out_ds.SetProjection(in_ds.GetProjection())
    geotransform = list(in_ds.GetGeoTransform())

    # Edit the geotransform so pixels are 10 times smaller than in original file in 100 m 
    geotransform[1] /= 10
    geotransform[5] /= 10
    out_ds.SetGeoTransform(geotransform)

    data = in_band.ReadAsArray(buf_xsize=out_columns, buf_ysize=out_rows)  # Specify a larger buffer size when reading data
    share = data / 10000 # Calculate the share of built-up per 100 m cell
    out_band = out_ds.GetRasterBand(1)
    out_band.SetNoDataValue(65535)
    out_band.WriteArray(share)

    out_band.FlushCache()
    out_band.ComputeStatistics(False)
    #out_ds.BuildOverviews('average', [2, 4, 8, 16, 32, 64])
    del out_ds
    del in_ds

# Function to read the template file geoinformation
def template_geoinfo(template_src):
    # Open the template and the input file 
    ras_ds_dst = gdal.Open(template_src)
    # Get the template geoinformation
    geot = ras_ds_dst.GetGeoTransform()
    xres_dst = geot[1]
    yres_dst = -geot[5]
    ULx_dst = geot[0]
    ULy_dst = geot[3]
    width_dst = int(ras_ds_dst.RasterXSize * xres_dst) # in meters
    height_dst = int(ras_ds_dst.RasterYSize * yres_dst) 
    # Get the projection
    crs_dst = ras_ds_dst.GetProjection()
    
    del ras_ds_dst
    
    return xres_dst, yres_dst, ULx_dst, ULy_dst, width_dst, height_dst, crs_dst

# Warp GHS data to 100 km tiles extent, with BUA resolution and projection
def warp_GHS(input_file, output_file, resample_alg, out_bounds, xres, yres, crs):
    warpOptions = gdal.WarpOptions(format = 'GTiff',
                                workingType = gdal.GDT_Float32,
                                outputType  = gdal.GDT_UInt16,
                                xRes = xres,
                                yRes = yres,
                                outputBounds  = out_bounds, #minX, minY, maxX, maxY
                                outputBoundsSRS  = crs,
                                dstSRS =crs,
                                dstNodata = 65535,
                                creationOptions = ['TFW=YES', 'COMPRESS=LZW', 'BIGTIFF=YES'],
                                resampleAlg = resample_alg)
                
    gdal.Warp(output_file, input_file, options=warpOptions)

# Function to run the preprocessing of GHS data - extraction and oversampling
# Takes the list of shapefiles as an input
# Return a list of corrupt tiles (shapefile paths)
def GHS_preprocessing(tiles_list): #list of shapefiles
    print('\nnumber of tiles:', str(len(tiles_list)))
    tiles_str = "\t".join(tiles_list)

    # 1. Save every tile as a seperate shapefile
    print('\nA. Save every tile as a seperate shapefile')
    tiles = gpd.read_file(tiles_USA_file)
    # Loop over each tile and save as a seperate shapefile
    for index, row in tiles.iterrows():
        key = (row.location).split(".")[0]
        key = key.split("/")[1]
        # If the tile is in the list of tiles to be processed, extract it and save
        if key+'/' in tiles_str: # check if the tile folder directory is in the list of tiles to be processed, hence +'/'
            print('Extract shapefile:', key)
            if not os.path.exists(os.path.join(root,'preprocessed-data',key)):
                os.makedirs(os.path.join(root,'preprocessed-data',key))
            path = os.path.join(root,'preprocessed-data',key,key+'.shp')
            if os.path.exists(path):
                os.remove(path) 
            tiles.loc[index:index].to_file(path)

    # 2. Oversample GHS data per tile
    print('\nB. Oversample GHS data per tile')
    # Assign tiles to be processed
    proc_list = tiles_list   
    # Extract GHS data to match the BUA extent
    for i, adate in enumerate([1975, 1990, 2000, 2015]):
        #For each tile, run the extraction
        for atile in proc_list: 
            aname = (os.path.basename(atile)).split(".")[0]
            print('Extract GHS in 1 m:', aname)
            try:
                # Extract GHS data for each tile
                filename = 'GHS'+str(adate)+'_'+aname+'_extract.tif'
                extract_file = os.path.join(root,'preprocessed-data',aname,filename) 
                extract_GHS(GHS_files[i], extract_file, atile)

                # Assign output path
                outfilename = 'GHS'+str(adate)+'_'+aname+'_10m.tif'
                oversample_file = os.path.join(root,'preprocessed-data',aname,outfilename)

                # Resample each tile to 10 m
                oversample_GHS(extract_file, oversample_file)

            except:
                print('CORRUPT:', aname)

# Function to check the results of processing an re-run the corrupt tiles
# Return the tiles still corupt after reprocessing
def check_reprocess(tiles_input_list):
    # Create an empty list to catch the corrupt tiles
    corrupt_list = []
    # Iterate the tiles per year and per tile
    for i, adate in enumerate([1975, 1990, 2000, 2015]):
        print(adate,': check tiles')
        for atile in tiles_input_list: 
            # Assign tif path
            aname = (os.path.basename(atile)).split(".")[0]
            outfilename = 'GHS'+str(adate)+'_'+aname+'_10m.tif'
            oversample_file = os.path.join(root,'preprocessed-data',aname,outfilename)
            # Check the input files
            # 1. Check for file size
            if os.stat(oversample_file).st_size == 0:
                corrupt_list.append(atile)
            elif os.stat(oversample_file.replace(".tif",".tfw")).st_size == 0:
                corrupt_list.append(atile)
            else:
                # 2. Try for gdal info
                try:
                    gdal_test = gdal.Info(oversample_file)
                except:
                    print('corrupt file', oversample_file)
                    corrupt_list.append(atile)

    corrupt_list = list(set(corrupt_list))
    # Re-process corrupt tiles, if any
    if len(corrupt_list)>0:
        print('\ncorrupt tiles:', corrupt_list)
        print('re-run corrupt tiles')
        GHS_preprocessing(corrupt_list)
        print('corrupt tiles re-processed')
    else:
        print('no corrupt tiles')
    
    return corrupt_list

################# DATA PATHS #################
# Load data
root = '/eos/jeodpp/data/projects/GHSL/2023_USA_MTVAL/'

# Load file with the AOI boundaries and uncertainity
county_file = os.path.join(root,'data','HISDAC-US','Cnty_Unc','County_Unc.shp')
#county_file = os.path.join(root,'data','HISDAC-US','Cnty_Unc','County_Unc_test.shp')

# GHS-MT layers
GHS_dir = os.path.join(root,'data','GHSL','R2023A')
GHS1975_file = os.path.join(GHS_dir, 'GHS_BUILT_S_E1975_GLOBE_R2023A_54009_100_V1_0.tif')
GHS1990_file = os.path.join(GHS_dir, 'GHS_BUILT_S_E1990_GLOBE_R2023A_54009_100_V1_0.tif')
GHS2000_file = os.path.join(GHS_dir, 'GHS_BUILT_S_E2000_GLOBE_R2023A_54009_100_V1_0.tif')
GHS2015_file = os.path.join(GHS_dir, 'GHS_BUILT_S_E2015_GLOBE_R2023A_54009_100_V1_0.tif')
GHS_files = [GHS1975_file, GHS1990_file, GHS2000_file, GHS2015_file]

# HISDAC-US layers
BUA_dir = os.path.join(root,'data','HISDAC-US','BUA')
BUA1975_file = os.path.join(BUA_dir, 'BUA_1975.tif')
BUA1990_file = os.path.join(BUA_dir, 'BUA_1990.tif')
BUA2000_file = os.path.join(BUA_dir, 'BUA_2000.tif')
BUA2015_file = os.path.join(BUA_dir, 'BUA_2015.tif')
BUA_files = [BUA1975_file, BUA1990_file, BUA2000_file, BUA2015_file]

OPPU_dir = os.path.join(root,'data','HISDAC-US','OPPU_1975_2015')
OPPU1975_file = os.path.join(OPPU_dir, 'data','OPPU_BUI_1975.tif')
OPPU1990_file = os.path.join(OPPU_dir, 'data','OPPU_BUI_1990.tif')
OPPU2000_file = os.path.join(OPPU_dir, 'data','OPPU_BUI_2000.tif')
OPPU2015_file = os.path.join(OPPU_dir, 'data','OPPU_BUI_2015.tif')
OPPU_files = [OPPU1975_file, OPPU1990_file, OPPU2000_file, OPPU2015_file]

# Ancillary data
GHS_mask_file = os.path.join(root,'data','GHSL','production-data','VRT','MT-production-mask.vrt')
tiles_file = os.path.join(root,'data','GHSL','tiles','tile_index.shp')

# Output working files
counties_MWD_file = os.path.join(root,'preprocessed-data','counties_MWD.shp')
tiles_USA_file = os.path.join(root,'preprocessed-data','tiles_USA.shp')

# Prepare directories
for f in ['preprocessed-data', 'processed-data']:
    if not os.path.exists(os.path.join(root,f)):
        os.mkdir(os.path.join(root,f))

if not os.path.exists(os.path.join(root,'preprocessed-data','vrt')):
    os.mkdir(os.path.join(root,'preprocessed-data','vrt'))

################# MAIN #################

# ### Preprocess USA counties
# print('\n1. Preprocess USA counties')
# print('- rasterize')
# # Rasterize counties and match to BUA 1975 file. The output will be the template for further processing.
# output_file = os.path.join(root,'processed-data','counties_raster.tif')
# rasterize_counties(county_file, output_file, BUA1975_file)
# # Use the county vector file as the boundary
# counties = gpd.read_file(county_file)
# # Extract boundary only
# print('- extract boundary')
# boundary = counties.dissolve()
# boundary.to_file(os.path.join(root,'processed-data','boundary.shp'))
# # Reproject counties to Mollweide projection, to be used for GHS preprocessing  
# print('- reproject to Mollweide')
# reproject_vector(county_file, 'ESRI:102039', counties_MWD_file, 'ESRI:54009')

### Select AOIs
# print('\n 2. Select the tiles and save as a list')
# # Read in GHS tiling schema in MWD resolution
# all_tiles = gpd.read_file(tiles_file)
# # Get the counties geometry in MWD
# counties_MWD = gpd.read_file(counties_MWD_file)
# # Dissolve
# print('- assign AOI')
# counties_MWD_diss = counties_MWD.dissolve()
# # Use the dissolved county file as the AOI boundary ans subset the tiles that cover AOI
# print('- subset tiles')
# tiles = all_tiles[all_tiles.intersects(counties_MWD_diss.geometry[0])]
# tiles.to_file(tiles_USA_file)
# ntiles = str(len(tiles))
# print('- number of tiles:', ntiles)
# # Loop over each tile and generate the path that will be used as the iterator
# tiles_list = []
# for index, row in tiles.iterrows():
#     key = (row.location).split(".")[0]
#     key = key.split("/")[1]
#     fileName = key+".shp"
#     if not os.path.exists(os.path.join(root,'preprocessed-data',key)):
#         os.makedirs(os.path.join(root,'preprocessed-data',key))
#     path = os.path.join(root,'preprocessed-data',key,fileName)
#     tiles_list.append(path)

# # Assign file to store tile_list
# tiles_csv = os.path.join(root,'preprocessed-data','tiles_list.csv')
# # Save list as dataframe
# tiles_df = pd.DataFrame(data={"tile_files": tiles_list})
# # Export to csv
# tiles_df.to_csv(tiles_csv, sep=',',index=False, header=False)

# # ############################## HERE IS THE sPACE FOR THE LOOP
# print('\n3. Oversample GHS data')
# Run GHS preprocessing: oversampling to 10 m resolution    
# GHS_preprocessing(tiles_list)


# ############################################################################
# Make this into a seperate script -->
# ############################################################################

print('\n3. Check the resampled GHS data for corrupt tiles')
# Read file storing the list of tiles
tiles_csv = os.path.join(root,'preprocessed-data','tiles_list.csv')
# Read as dataframe
tiles_df = pd.read_csv(tiles_csv, sep=',', header=None)
# Flatten into a list
tiles_list = [ item for sublist in tiles_df.values.tolist() for item in sublist ]
#tiles_list = tiles_list[0:50]
print('Number of tiles:',str(len(tiles_list)),'\n')
# Check the resampled tiles, re-run if necessery
check_list=tiles_list
while len(check_list)>9: # tiles less than 10 can be handled manually
    print('check for corrupt tiles:', str(len(check_list)),'tiles')
    check_list = check_reprocess(check_list)

print('\n4. Reproject GHS data')
# Select the rasterized counties as the template for the grid size, position and projection
template_file = os.path.join(root,'processed-data','counties_raster.tif')
# Get template geoinfo
xres_dst, yres_dst, ULx_dst, ULy_dst, width_dst, height_dst, crs_dst = template_geoinfo(template_file)
# Extract GHS data to match the BUA extent
for i, adate in enumerate([1975, 1990, 2000, 2015]):
    print(adate)
    # Make a list for vrt paths
    GHS_vrt_files = []
    #For each tile, run the extraction
    print('Collect tiff paths...')
    for atile in tiles_list: 
        # Assign tif path
        aname = (os.path.basename(atile)).split(".")[0]
        outfilename = 'GHS'+str(adate)+'_'+aname+'_10m.tif'
        oversample_file = os.path.join(root,'preprocessed-data',aname,outfilename)
        GHS_vrt_files.append(oversample_file)

    # Build vrt
    print('Build vrt', str(adate))
    vrt_file = os.path.join(root,'preprocessed-data','vrt','GHS'+str(adate)+'_10m.vrt')
    gdal.BuildVRT(vrt_file, GHS_vrt_files, options=gdal.BuildVRTOptions(outputSRS='ESRI:54009', 
                                                                        srcNodata=65535,
                                                                        VRTNodata=65535))

    # Now save the vert file as the tif file, matched to BUA extent and resolution
    print('Build final tiff file', str(adate))
    outbounds_dst = [ULx_dst,ULy_dst-height_dst,(ULx_dst+width_dst),(ULy_dst)]
    output_file = os.path.join(root,'processed-data','GHS'+str(adate)+'_warped.tif') 
    warp_GHS(vrt_file, output_file, 'sum',outbounds_dst, xres_dst, yres_dst,crs_dst)

print('done')


                
###############     NOT USED    ######################

# # Create a 100 km tiling schema
# dx = 100 * 1000 #n*[1km]
# dy = 100 * 1000 #n*[1km]
# nrows = math.ceil(height_dst / dy)
# ncols = math.ceil(width_dst / dx)
# print('rows,cols:',str(nrows),str(ncols))

# # Make a list for vrt paths
# GHS_vrt_files = []
# print('Warp input data to match the 100 km block and BUA resolution')
# for i, adate in enumerate([1990]):#, 1990, 2000, 2015]):

#     # Get vrt with 10m GHS data
#     vrt_name = 'GHS'+str(adate)+'_10m.vrt'
#     vrt_file = os.path.join(root,'preprocessed-data',vrt_name)
#     GHS_vrt_files.append(vrt_file)
#     # Warp GHS data
#     output_file_base = os.path.join(root,'preprocessed-data','warped','GHS'+str(adate)) 
#     # Create a set of 100 km tiles
#     i=0
#     warped_tiles = []
#     corrupt_warps = []
#     for nrow in range(nrows):
#         for ncol in range(ncols):
#             try:
#                 print(str(i), str(nrow),str(ncol))
#                 ULx = ULx_dst + dx*ncol
#                 ULy = ULy_dst - dy*nrow
#                 outbounds = [ULx,ULy-dx, (ULx+dy),(ULy)],
#                 output_file = output_file_base+'_'+str(i)+'.tif'
#                 warped_tiles.append(output_file)
#                 warp_GHS(vrt_file, output_file, 'sum', outbounds, xres_dst, yres_dst,crs_dst)
#                 i=i+1
#             except:
#                 print('fail:', output_file)

#             # Check the results
#             try:
#                 gdal_test = gdal.Info(output_file)
#             except:
#                 print('warp fail', output_file)
#                 corrupt_warps.append([nrow, ncol,output_file])
    
#     # Build a vrt
#     vrt_file = output_file_base+'_warped.vrt'
#     gdal.BuildVRT(vrt_file, warped_tiles, options=gdal.BuildVRTOptions(outputSRS=crs_dst, 
#                                                                          srcNodata=65535,
#                                                                          VRTNodata=65535))
    
#     print(str(adate),'GHS vrt built') 

#     # Now save the vert file as the tif file, matched to BUA extent and resolution
#     outbounds_dst = [ULx_dst,ULy_dst-height_dst,(ULx_dst+width_dst),(ULy_dst)]
#     output_file = os.path.join(root,'processed-data','GHS'+str(adate)+'_match.tif') 
#     warp_GHS(vrt_file, output_file, 'near',outbounds_dst, xres_dst, yres_dst,crs_dst)


# ovsamlpef = '/eos/jeodpp/data/projects/GHSL/2023_USA_MTVAL/preprocessed-data/tile_r39_c108_6/GHS1975_tile_r39_c108_6_10m.tif'    
# gdal_test = gdal.Info(ovsamlpef)
# print(gdal_test)
