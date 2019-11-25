import pandas as pd

#PURPOSE: calculate reach averaging for widths not tagged with lakes/resevoirs, and distance <100m from its closest MERIT reach
#INPUT: all table.csv calculated from Step 1
#OUTPUT: a csv file containing reach-average width values, and count (count <30 are discarded for machine learning training)

df = pd.read_csv('all_table_MERIT_GRWL.csv')

print(len(df))
df = df[df['lakeFlag']==0]
df = df[df['distance']<=0.001]
print(len(df))

data = df.groupby('COMID')['width_m'].mean().reset_index()  #averaging by COMID
data['count']=df.groupby('COMID')['distance'].count().values
data['elev_m'] = df.groupby('COMID')['elev_m'].mean().reset_index()['elev_m'].values
data.to_csv('reach_average_GRWL.csv',index=False)
