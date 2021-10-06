#!/usr/bin/env python
# coding: utf-8

# In[9]:


import pandas as pd
import geopandas as gpd
from shapely.geometry import Polygon, LineString, Point
from numpy import random
from osgeo import gdal, ogr


#def folder paths
#still need to change pathfile in function itself for oblast selection
def folder_paths():
    print('Folder Paths Function')
    default = input('Default folder_paths? (Yes or no): ')
    if default == 'Yes':
        #set path for original raster masks
        og_rasters_folder = '/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/Mask/Russia/'
        #set path for storing ediited rasters
        raster_process_folder = '/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/Mask/vectorized_mask/raster_process/'
        #set path for storing polygonized data and polygons
        poly_process_folder = '/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/Mask/vectorized_mask/polygonized/'
        #final outputs
        final_outputs = '/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/Mask/vectorized_mask/final_outputs/'
    else: 
        og_rasters_folder = input('Enter filepath for folder containing original masks:') + '/'
        raster_process_folder = input('Enter filepath for empty folder to store processed rasters ') + '/'
        poly_process_folder = input('Enter filepath of empty folder to store processed shapefiles') + '/'
        final_outputs = input('Enter filepath of empty folder to store final outputs of data') + '/'
    #add .tif for saving rasters
    tif_save = '.tif'
    #add .shp for saving shapefiles
    shp_save = '.shp'
    print('Folder Paths Completed')
    print('|----------------------------------------|')
    return og_rasters_folder, raster_process_folder,poly_process_folder,final_outputs, tif_save, shp_save

#Select parameters function
def select_param(tif_save):
    print('Select Parameters Function')
    years_ls = []
    years_in = ''
    mask_name = 'Russia_'
    default_param = input('Use default parameters for model? (Yes or no)')
    if default_param == 'Yes':
        id_select = 77
        filter_size_init = 6000
        filter_size_overlay = 250
        samples = 10
        years_ls = ['Russia_2020.tif','Russia_2019.tif','Russia_2018.tif','Russia_2017.tif','Russia_2016.tif']
    else:
        mask_name = input('Enter mask name') + '_'
        print('Enter "stop" to continue')
        while years_in != 'stop':
            years_in = input('Enter year of mask to process: ')
            if years_in != 'stop':
                years_in = mask_name + years_in + tif_save
                years_ls.append(years_in)
        ##select oblast/region
        id_select = int(input('Enter ID of which oblast to select: '))
        ##select filters for polygon area
        print('Enter default or enter integer of area in m^2/1000 to filter size of polygons for initial and overlay polygons')
        filter_size_init= input('Enter size of poly_mask filter: ')
        #select number of random point samples per year
        samples = input('Enter number of samples to collect per polygon year: ')
        if filter_size_init == 'default':
            filter_size_init = 6000
        else:
            filter_size_init = int(filter_size_init)
        filter_size_overlay = input('Enter size of overlay filter: ')
        if filter_size_overlay == 'default':
            filter_size_overlay = 250
        else:
            filter_size_overlay = int(filter_size_overlay)
        print('Parameter Selection Completed')
        print('|----------------------------------------|')
    return years_ls, id_select, samples, filter_size_init, filter_size_overlay

#def select oblast to analyze
def oblast_select(id_select):
    print('Oblast Select Function')
    print('ID Options: 77 = Volgograd ', )
    ukr_rus_dat = pd.read_csv("/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/master/ukr_rus_master.csv")
    rus_shp = gpd.read_file("/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/data/shapefile/Russia_admin_clean/rus_admin_cl.shp")
    ukr_shp = gpd.read_file("/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/data/shapefile/Ukraine_admin/regions.shp")
    ukr_shp['ID'] = ukr_shp['ID'] + 100
    ukr_shp = ukr_shp[['ID','geometry']]
    rus_shp = rus_shp.rename(columns = {"ID_1":"ID"})
    rus_shp = rus_shp[['ID','geometry']]
    shp_stack = [rus_shp,ukr_shp]
    shp = pd.concat(shp_stack)
    #merge and set geometry of master dataset
    x = ukr_rus_dat.merge(shp, on = 'ID').set_geometry('geometry')
    m = x[['ID','Country','Oblast-eng','geometry']].drop_duplicates()
    #feature collection for ee conversion
    selected_oblast = m[m['ID'] == id_select]
    selected_oblast.to_file('/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/Mask/vectorized_mask/shapefile_select/selected_oblast.shp')
    print('Oblast Select Completed')
    print('|----------------------------------------|')
    return selected_oblast


#function for reprojecting/warping original mask
def warp_raster(years_ls,og_rasters_folder,raster_process_folder,tif_save):
    from osgeo import gdal
    print('Warp Raster function')
    print('Set to convert to 4326 projection')
    for i in years_ls:
        #enter directory folder
        folder_input  = og_rasters_folder
        folder_output = raster_process_folder
        raster_name   =  i
        #remove .tif
        i = i.replace('.tif','')
        output_raster = i + '_warp'
        mask_input    = folder_input + raster_name 
        mask_output   = folder_output + output_raster + tif_save
        #open raster
        input_raster  = gdal.Open(mask_input)
        warp = gdal.Warp(mask_output,mask_input,dstSRS='EPSG:3857')
    print('Warp Raster Completed')
    print('|----------------------------------------|')
    return warp

#function for clipping winter wheat mask to shapefile
def clip_raster_with_shape(years_ls,raster_process_folder,tif_save):
    print('Clip raster with shape function')
    from osgeo import gdal, ogr
    for i in years_ls:
        folder_input = raster_process_folder
        folder_output = raster_process_folder
        raster_name =  i.replace('.tif','') + '_warp.tif'
        input_raster = gdal.Open(folder_input + raster_name)
        output_name = raster_name.replace('_warp.tif','') + '_clipped' + tif_save
        output_raster = folder_output + output_name
        shape_clip = '/Users/christianabys/Desktop/School/Maryland/Research/Wheat_Forecast/Mask/subselect_shapefile/oblast_select.shp'
        OutTile = gdal.Warp(output_raster,input_raster,cutlineDSName = shape_clip,dstSRS='EPSG:4326')
    print('Clip Raster Completed')
    print('|----------------------------------------|')
    return OutTile


#EROSION FUNCTION TO DO 

#function for polygonizing clipped winter wheat mask
#function to polygonize mask, returns a shapefile projected to 4326
def polygonize_gdal(years_ls,raster_process_folder,tif_save,shp_save,poly_process_folder):
    from osgeo import gdal
    print('Polygonize raster with GDAL function')
    for i in years_ls:
        folder_input  = raster_process_folder
        raster_name   = i.replace('.tif','') + '_clipped' + tif_save
        input_raster  = folder_input  + raster_name
        input_raster  = gdal.Open(input_raster)
        folder_output = poly_process_folder
        folder_output_name = raster_name.replace('_clipped.tif','') + '_polygonized' + shp_save
        output_shp = folder_output + folder_output_name
        drv = ogr.GetDriverByName('ESRI Shapefile')
        outfile = drv.CreateDataSource(output_shp)
        band = input_raster.GetRasterBand(1)
        outlayer = outfile.CreateLayer('polygonized raster', srs = None )
        newField = ogr.FieldDefn('DN', ogr.OFTReal)
        outlayer.CreateField(newField)
        gdal.Polygonize(band, None, outlayer, 0, [])
        outfile = None
        shapefile = gpd.read_file(output_shp)
        shapefile = shapefile.set_crs('EPSG:4326')
        shapefile.to_file(output_shp)
    print('Polygonize Raster Completed')
    print('|----------------------------------------|')
    return shapefile


#define function for filtering out unneeded DN,calc area, and divide by 1000
#polygonized mask, year of mask (TO DO: code in grabbing year of mask)
def mask_poly_calc(years_ls,poly_process_folder,shp_save,filter_size_init):
    print('Calculate Polygon Mask Function')
    for i in years_ls:
        folder_input = poly_process_folder
        input_name = i.replace('.tif','') + '_polygonized' + shp_save
        output_name = i.replace('.tif','') + '_poly_filt' + shp_save
        poly_mask = gpd.read_file(folder_input + input_name)
        poly_mask = poly_mask.to_crs("EPSG:3857")
        poly_mask = poly_mask[poly_mask['DN'] != 0]
        poly_mask['m2/1000'] = poly_mask.area/1000
        poly_mask = poly_mask.drop(columns = ['DN'])
        poly_mask = poly_mask[poly_mask['m2/1000'] > int(filter_size_init)]
        poly_mask['year'] = int(i.replace('Russia_','').replace('.tif',''))
        #create column for storing the year that the mask represents
        poly_mask = poly_mask.to_crs('EPSG:4326')
        poly_mask.to_file(folder_input + output_name) 
    print('Calculate Polygon Mask Completed')
    print('|----------------------------------------|')


#function for overlaying two masks over one another and getting areas NOT intersected
##additonally reorganizes columns
##explodes polygons into individual units
##recalculates exploded polygons
def overlay_mask(mask1,mask2,filter_size_overlay):
    print('Overlay Mask: Symmetric Differences')
    #conduct overlay
    mask = gpd.overlay(mask1,mask2, how = 'symmetric_difference')
    mask = mask.to_crs('EPSG:3857')
    #fillna in columns
    mask.iloc[:, [0,1,2,3]] = mask.iloc[:, [0,1,2,3]].fillna(0)
    #combine columns
    mask['year'] = mask['year_1'].map(int) + mask['year_2'].map(int)
    mask['m2/1000'] = mask['m2/1000_1'].map(float) + mask['m2/1000_2'].map(float)
    #drop unneeded columns
    mask = mask.iloc[:,[4,5,6]]
    #explode geometries
    mask = mask.explode()
    #recalculate area
    mask['m2/1000'] = mask.area/1000
    #filter out small areas
    mask = mask[mask['m2/1000'] > int(filter_size_overlay)]
    #reindex
    mask = mask.reset_index()
    #drop columns
    mask = mask.iloc[:,[2,3,4]]
    mask = mask.to_crs('EPSG:4326')
    print('Overlay Mask Completed')
    print('|----------------------------------------|')
    return mask

#function to read in masks and apply to overlay_mask
def mask_create(years_ls,poly_process_folder,shp_save):
    print('Create Mask Function')
    #sort list ascending
    years_ls.sort(reverse=True)
    #seperate out mask files
    mask1,mask2,mask3,mask4,mask5 = years_ls
    #retrieve individual mask shapefiles
    mask1 = mask1.replace('.tif','_poly_filt') + shp_save
    mask1 = gpd.read_file(poly_process_folder + mask1)
    mask2 = mask2.replace('.tif','_poly_filt') + shp_save
    mask2 = gpd.read_file(poly_process_folder + mask2)
    mask3 = mask3.replace('.tif','_poly_filt') + shp_save
    mask3 = gpd.read_file(poly_process_folder + mask3)
    mask4 = mask4.replace('.tif','_poly_filt') + shp_save
    mask4 = gpd.read_file(poly_process_folder + mask4)
    mask5 = mask5.replace('.tif','_poly_filt') + shp_save
    mask5 = gpd.read_file(poly_process_folder + mask5)
    print('Create Mask Function Completed')
    print('|----------------------------------------|')
    return mask1,mask2,mask3,mask4,mask5

#run mask overlay 
def run_mask(mask1,mask2,mask3,mask4,mask5,filter_size_overlay,final_outputs,shp_save):
    print('Run Mask Function')
    mask6 = overlay_mask(mask1,mask2,filter_size_overlay)
    mask7 = overlay_mask(mask3,mask4,filter_size_overlay)
    mask8 = overlay_mask(mask6,mask7,filter_size_overlay)
    mask_final = overlay_mask(mask8,mask5,filter_size_overlay)
    mask_final.to_file(final_outputs + 'mask_final' + shp_save)
    print('Run Mask Completed')
    print('|----------------------------------------|')
    return mask_final

#select random points in each polygon at n samples and return point with associated year
def random_points(mask_final,samples,shp_save,final_outputs):
    print('Random Points in Polygon')
    import time
    #dissolve all polygons by year and reset index
    mask = mask_final.dissolve(by = 'year').reset_index()
    #initalize point and year ls
    points = []
    year_ls = []
    #run through each row
    for index, row in mask.iterrows():
            #grab year of polygon
            year = row[0]
            #initalize counter
            i = 0
            sub_mask = mask[mask['year'] == year]
            min_x, min_y, max_x, max_y = sub_mask.bounds.iloc[0]
            print('Subset Mask Year: ',sub_mask.year)
            while i < samples:
                start_time = time.time()
                point = Point(random.uniform(min_x, max_x), random.uniform(min_y, max_y))
                if sub_mask.contains(point).bool():
                    #append points to point list
                    points.append(point)
                    #append year to year list
                    year_ls.append(year)
                    #add to counter
                    i += 1
                    print('Sample: ', i)
                    #add points,year to a dataframe
                    point_df = pd.DataFrame(list(zip(points,year_ls)))
                    print(point_df)
                    print("Execution Time:", "%s seconds" % (time.time()-start_time))
                    print('|--------------------------------------------------------|')

    #convert points from list to geodataframe and assign geometry and crs
    point_df = point_df.rename(columns = {0:'geometry',1:'year'})
    point_gdf = gpd.GeoDataFrame(point_df,geometry = 'geometry',crs = "EPSG:4326")
    point_gdf.to_file(final_outputs + 'sampled_points' + shp_save)
    print('Random Points in Polygon Sampling Completed')
    print('|--------------------Sampling-Routine-Completed--------------------|')
    return point_gdf

#define routine
def main():
    og_rasters_folder, raster_process_folder,poly_process_folder,final_outputs, tif_save, shp_save = folder_paths()
    years_ls, id_select, samples, filter_size_init, filter_size_overlay = select_param(tif_save)
    oblast_select(id_select)
    warp_raster(years_ls,og_rasters_folder,raster_process_folder,tif_save)
    clip_raster_with_shape(years_ls,raster_process_folder,tif_save)
    polygonize_gdal(years_ls,raster_process_folder,tif_save,shp_save,poly_process_folder)
    mask_poly_calc(years_ls,poly_process_folder,shp_save,filter_size_init)
    mask1,mask2,mask3,mask4,mask5 = mask_create(years_ls,poly_process_folder,shp_save)
    mask_final = run_mask(mask1,mask2,mask3,mask4,mask5,filter_size_overlay,final_outputs,shp_save)
    random_points(mask_final,samples,shp_save,final_outputs)
#run main routine
if __name__ == "__main__":
        main()
#function to grab pixel values from each point
#tensor flow function for image classification
#function to compute statistics from each inidivudal pixel value


# In[ ]:




