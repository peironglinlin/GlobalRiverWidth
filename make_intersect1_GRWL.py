import glob
from osgeo import ogr
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point,Polygon

#intersect MERIT width GRWL such that reach-averaging in Step 2 can be conducted
#INPUT: all MERIT flowlines shapefile; all GRWL shapefiles
#OUTPUT: a csv file containing their join [segmentID,segmentInd,COMID,distance,width_m,nchannels,lakeFlag,lon,lat,elev_m]

def generate_MERIT_bounds(files):
    names = []
    latS = []
    latN = []
    lonW = []
    lonE = []
    for i in range(0,nf):
        print('... reading '+files[i]+' ...')
    #     data = gpd.read_file(files[i])  #do not use geopandas
    #     import pdb;pdb.set_trace()
        ds = ogr.Open(files[i])
        layername = files[i].split('.')[0].split('/')[-1]
        print(layername)
        lyr = ds.GetLayerByName(layername)
        lyr.ResetReading()

        latss = []
        latnn = []
        lonww = []
        lonee = []
        for feat in lyr:
            env = feat.GetGeometryRef().GetEnvelope() # minx, miny, maxx, maxy format
            latss.append(env[2])
            latnn.append(env[3])
            lonww.append(env[0])
            lonee.append(env[1])

        #compute bounding box
        latS.append(min(latss))
        latN.append(max(latnn))
        lonW.append(min(lonww))
        lonE.append(max(lonee))
        names.append(layername)
        print((latS[-1],latN[-1],lonW[-1],lonE[-1]))
    
    #create final dataframe to write to file
    df = pd.DataFrame({'name':names,'latS':latS,'latN':latN,'lonW':lonW,'lonE':lonE})
    return df

def make_polygons(df):
    npoly = len(df)
    polygons = []
    for i in range(0,npoly):
        lonlist = [df.lonW[i],df.lonE[i],df.lonE[i],df.lonW[i],df.lonW[i]]
        latlist = [df.latS[i],df.latS[i],df.latN[i],df.latN[i],df.latS[i]]
        polygons.append(Polygon(zip(lonlist,latlist)))
    return polygons

def join_MERIT_with_GRWL(files, polygons):
    grwl = []
    regions = []
    for i in range(0,nf):  #loop through all GRWL vectors
        print('*********************************')
        print(i)
        data = gpd.read_file(files[i])
        now = files[i].split('/')[-1].split('.')[0]
        print(now)
    #     if now == 'SE32':
        nfeature = len(data)
        if nfeature != 0:
            cc = gpd.GeoDataFrame(geometry=data.geometry.centroid)  #centerline of measurements (every 30-50m)
            lon = cc.geometry.x
            lat = cc.geometry.y
            latN = max(lat)
            latS = min(lat)
            lonW = min(lon)
            lonE = max(lon)
            lonll = [lonW,lonE,lonE,lonW,lonW]
            latll = [latS,latS,latN,latN,latS]
            poly = Polygon(zip(lonll,latll))

            for j in range(0,npoly):
                if poly.intersects(polygons[j]):
                    print(df.name[j])
                    grwl.append(now)
                    regions.append(df.name[j])
            print('*********************************')
        else:
            grwl.append(now)
            regions.append('-999')
    dfnew = pd.DataFrame({'GRWL':grwl,'MERIT':regions})
    return dfnew

def final_join_MERIT_GRWL(df):
    nf = len(df)
    allpoints = pd.DataFrame({})
    for i in range(0,nf):
        if df.MERIT[i] != '-999':
            fgrwl = '../data/GRWL_vector_V01.01/'+df.GRWL[i]+'.shp'
            fmerit = '/tigress/peirongl/MERIT/raster/cleaned/level_01/'+df.MERIT[i]+'.shp'
            print('********************************')
            print('   GRWL = '+fgrwl)
            print('   MERIT = '+fmerit)
            meritnow = gpd.read_file(fmerit)

            #read GRWL and create centroid point
            grwl = gpd.read_file(fgrwl)
            nfeature = len(grwl)
            if nfeature != 0:
                cc = gpd.GeoDataFrame(grwl,geometry=grwl.geometry.centroid)  #centerline of measurements (every 30-50m)

                #create buffer
                print(' create buffer... wait ...')
                buffersize = 0.01 #~1km
                poly = cc.buffer(buffersize)
                polygpd = gpd.GeoDataFrame(cc[['width_m', 'nchannels', 'segmentID',\
                                               'segmentInd', 'lakeFlag', 'lon', 'lat', 'elev_m',]],geometry=poly)

                #spatial join
                print(' spatial join with flowlines.. wait ...')
                join = gpd.sjoin(polygpd,meritnow,how='inner',op='intersects')
                merge=join.merge(meritnow,on='COMID',how='left')
                print(' calculating distance.. wait ...')
                merge['distance']=[Point(merge['lon'][i],merge['lat'][i]).distance(merge['geometry_y'][i]) for i in range(0,len(merge))] 
                join11 = merge.groupby(['segmentID','segmentInd']).agg({'distance':'min'}).reset_index() #min dist: width and MERIT
                merge11 = join11.merge(merge,on=['segmentID','segmentInd','distance'],how='left')
        #         import pdb;pdb.set_trace()
                final = merge11[['segmentID','segmentInd','COMID','distance','width_m','nchannels','lakeFlag','lon','lat','elev_m']]
                final.to_csv('tmp/'+df.GRWL[i]+'_'+df.MERIT[i]+'.csv',index=False)
                allpoints = allpoints.append(final)
    return allpoints

#list all MERIT flowline shapefiles
files = sorted(glob.glob('/tigress/peirongl/MERIT/raster/cleaned/level_01/*riv*.shp'))
nf = len(files)
df = generate_MERIT_bounds(files)
polygons = make_polygons(df)

#list all GRWL river width
files = glob.glob('../data/GRWL_vector_V01.01/*.shp')
nf = len(files)
dfnew = join_MERIT_with_GRWL(files, polygons)

#final join GRWL width MERIT
allpoints = final_join_MERIT_GRWL(dfnew)
fon = 'all_table_MERIT_GRWL.csv'
print('... writing to '+fon+' ...')
allpoints.to_csv(fon,index=False)
