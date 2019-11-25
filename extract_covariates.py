import rasterio as rio
import pandas as pd
import numpy as np

#Input #1: a list of GeoTiff datasets containing the covariate values
#Input #2: GeoTiff datasets for global river unit catchments (polygon to raster, same resolution as Input #1)
#Output: a csv file containing zone ID (e.g., COMID here) and their covariate values

zones = 'COMID'
res = '001'

#data path for all covariates (pre-processed GeoTiff dataset, bilinearly interpolated to ~1km resolution globally)
all_data = [
    ('AI','../../RESEARCH_DATA/ai_et0/AI_0125_%s.tif'%res),  #aridity
    ('LAI','../../RESEARCH_DATA/lai_avhrr/LAI_0125_001.tif'), #LAI
    ('SND','../../RESEARCH_DATA/soilgrid250m/SNDPPT_001.tif'), #sand
    ('CLY','../../RESEARCH_DATA/soilgrid250m/CLYPPT_001.tif'), #clay
    ('SLT','../../RESEARCH_DATA/soilgrid250m/SLTPPT_001.tif'), #silt
    ('Dom','../../RESEARCH_DATA/wateruse/domww_001_20102014_mean.tif'), #domestic water use
    ('Ind','../../RESEARCH_DATA/wateruse/indww_001_20102014_mean.tif'), #industrial water use
    ('Irr','../../RESEARCH_DATA/wateruse/irrww_001_20102014_mean.tif'), #irrigational water use
    ('Urb','../../RESEARCH_DATA/urban_landsat/raster_urban2010.tif'), #urban fraction
    ('WTD','../../RESEARCH_DATA/WTD_Fan2013_science/WTD_0125_001.tif')] #water table depth
n_covar = len(all_data)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx #return the index of the nearest value

def mycount(array):
    return len(array[np.isnan(array)])/len(array)

# final = pd.DataFrame({zones:np.unique(catw)})
for ivar in range(n_covar):
    varnow = all_data[ivar][0]
    fin = all_data[ivar][1]
    print('   ... processing %s ...'%varnow)
    ds0 = rio.open(fin)
    #read data to be averaged
    nlat,nlon = ds0.shape
    xres = ds0.transform[0]
    yres = ds0.transform[4]
    lonww = ds0.transform[2]
    latnn = ds0.transform[5]
    lats = np.linspace(latnn,latnn+yres*nlat,nlat)
    lons = np.linspace(lonww,lonww+xres*nlon,nlon)
    
    output = pd.DataFrame({})
    #read data of zones
    for pfaf in range(1,9):
        print('      ... pfaf = %02d ...'%pfaf)
        ds = rio.open('../../RESEARCH_DATA/tif_merit/cat_pfaf_%02d.tif'%pfaf)  #input data: pre-processed unit catchment GeoTiff files (polygon to raster, ~1km resolution)
        lonW,latS,lonE,latN = ds.bounds
        nlat,nlon = ds.shape
        catw = ds.read(1)
        
        # find the indices of the latitude lower and upper index
        latmin = find_nearest(lats,latN)  #note the lat starts from above, so need to use latN
        lonmin = find_nearest(lons,lonW)
        latmax = latmin+nlat
        lonmax = lonmin+nlon

        #subsetting global data using the averaging zone's bounding box
        if varnow == 'AI':
            mydata = (ds0.read(1)/10000).astype('float')
        else:
            mydata = ds0.read(1).astype('float')
        #subsetting
        if pfaf != 3:
            cutdata = mydata[latmin:latmax , lonmin:lonmax]  #note lats start from 90N; only read subset
        else:
            cutdata1 = mydata[latmin:latmax , lonmin:lonmax]
            fixlon = nlon-cutdata1.shape[1]
            cutdata2 = mydata[latmin:latmax , 0:fixlon]
            cutdata = np.append(cutdata1,cutdata2,axis=1)
            del cutdata1,cutdata2
        print(cutdata.shape)

        #preprocessing covariates - unify missing values
        if (varnow == 'AI') | (varnow == 'WTD') | (varnow == 'topo'):
            mydata[mydata<0] = np.nan
        elif (varnow == 'litho'):
            mydata[mydata==0] = np.nan
        elif (varnow == 'LAI'):
            mydata[mydata<0] = 0. #LAI set as zero in nodata region
        elif (varnow == 'SND') | (varnow == 'CLY') | (varnow == 'SLT'):
            mydata[mydata>100] = np.nan
        else:
            print('')
        del mydata    

        #create dataframe
        dd = pd.DataFrame({zones:catw.flat,varnow:cutdata.flat})
        dd = dd[dd[zones]>0] #only average for places with COMID
        output = output.append(dd.groupby(zones)[varnow].mean().reset_index())  #if nan value exists, retain NaN
    
#     import pdb;pdb.set_trace()
    if ivar == 0:
        final = output
    else:
        final = final.merge(output,on=zones,how='left')
    del cutdata

fon='data/features_%s.csv'%(res)
print('... writing to '+fon+ ' ...')
final.to_csv(fon,index=False)
