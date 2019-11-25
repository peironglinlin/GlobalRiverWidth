import pandas as pd

#PURPOSE: conduct reach-averaging with >=30 valid readings
#INPUT: MERIT_width_table.csv calculated from Step 1
#OUTPUT: A csv file containing reach-averaged width

fin = 'MERIT_width_table.csv'
df = pd.read_csv(fin)

comids = []
www = []
nf = len(df)
for i in range(0,nf):
    print('... processing COMID = %s ...'%df.COMID[i])
    widths = [float(x) for x in df.elev[i].split('[')[1].split(']')[0].split(',')]
    ww = 0
    n = 0
    for j in widths:
        if j!=-9999:  #!!!only trust MERIT Hydro width when it is larger than
            n = n+1
            ww = ww+j
#     import pdb;pdb.set_trace()    
    if n>=30:
        ww = ww/n
        comids.append(df.COMID[i])
        print('   ww = %.2f'%ww)
        www.append(ww)

df0 = pd.DataFrame({'COMID':comids,'MERIT_elev':www})
fon = 'MERIT_width_table_averaged.csv'
print('... writing to %s ...'%fon)
df0.to_csv(fon,index=False)
    
