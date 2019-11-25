import geopandas as gpd
import pandas as pd
import numpy as np
import rasterio as rio
import math

#INPUT: 8 shapefiles containing the MERIT river flowlines and their intersected MERIT tile number (pre-processed using Python GeoPandas sjoin function)
#INPUT: all MERIT Hydro tiles (1150 tiles) stored in GeoTiff format
#OUTPUT: A csv file listing COMID, and all cross-sectional widths read from different MERIT Hydro tiles

#read spatial join table
df = gpd.GeoDataFrame({})
for pfaf in range(1,9):
    fin = 'intersect_tables/pfaf_%02d_merit_tile.shp'%pfaf
    df = df.append(gpd.read_file(fin))
df = df.dropna()
nf = len(df)

#how many number of MERIT tiles are there
nf0 = len(np.unique(df.tif))

res = 0.000833333333333
iii = 0
comids = []
mytif = []
mywidth = []
for i in range(0,nf0):  #only open nf0 times the tif file
    tilenow = np.unique(df.tif)[i].split('_')[0]
    tifnow = '../../../MERIT/width/tif/%s_wth.tif'%tilenow
#     if np.unique(df.tif)[i] == 's05w045_wth': #TESTING WITH s05w045 width file
    print('... openning %s ...'%tifnow)
    #open tifdata
    ds = rio.open(tifnow)
    lonW,latS,lonE,latN = ds.bounds
    data = ds.read(1)

    #read flowline coordinates and extract width value
    dfnow = df[df.tif == np.unique(df.tif)[i]]
    nfnow = len(dfnow)
    for j in range(0,nfnow):
#         import pdb;pdb.set_trace()
        print('    ... read COMID = %s ...'%dfnow['COMID'].iloc[j])
        comids.append(dfnow['COMID'].iloc[j])
        mytif.append(tifnow)
        widths = []
        #read polyline data
        coords = dfnow.geometry.iloc[j].coords
        for k in range(0,len(coords)):
#                 import pdb;pdb.set_trace()
            lonnow = coords[k][0]
            latnow = coords[k][1]
            indlon = math.floor(int((lonnow-lonW)/res))
            indlat = math.floor(int((latN-latnow)/res))
            try:
                ww = data[indlat,indlon]
            except:
                ww = -9999.
            print('   ... ww = %s ...'%ww)
            widths.append(ww)
        mywidth.append(widths)

df = pd.DataFrame({'COMID':comids,'tif':mytif,'width':mywidth})
fon = 'MERIT_width_table.csv'
print('... writing to %s ...'%fon)
df.to_csv(fon,index=False)
