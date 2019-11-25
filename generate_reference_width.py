import pandas as pd

#PURPOSE: calculate reference width by averaging GRWL and MERIT Hydro reach-averaged width
#INPUT: reach-averaged GRWL and reach-average MERIT Hydro width
#OUTPUT: a csv file containing their mean values at overlapping reaches

n=30

df1 = pd.read_csv('MERIT_width_table_averaged.csv')
df2 = pd.read_csv('reach_average_GRWL.csv')
df2 = df2[df2['count']>=n]

df = df1.merge(df2,on='COMID',how='inner').drop_duplicates()
print(len(df))

df['width_mean'] = (df['width_m']+df['MERIT_width'])/2
df = df[['COMID','width_m','MERIT_width','width_mean']]
df.columns = ['COMID','GRWL','MERIT','width_mean']
df = df.groupby('COMID')[df.columns.drop('COMID')].mean().reset_index()
df = df.sort_values(by='COMID')
df.to_csv('new_mean_width_n_%s.csv'%n,index=False)
