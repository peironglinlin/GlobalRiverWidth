from netCDF4 import Dataset
import pandas as pd
import numpy as np
import xarray as xr
import os
from scipy.stats import norm,lognorm,gumbel_r,pearson3
import scipy.stats
# import multiprocessing as mp
from mpi4py import MPI

#MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

ndays = 12784
time = pd.date_range(start='1/1/1979', periods=ndays, freq='D')
indices = [np.where(time=='%s-01-01'%yyyy)[0][0] for yyyy in range(1979,2014)]
distributions = [norm,lognorm,gumbel_r,pearson3] #normal, log-normal,gumbel,log-pearson Type III

def get_annual_max_q(qq):
    """return N years of annual maximum Q; input: 35-year of daily streamflow""" 
    qmax = [max(qq[index:indices[i+1]]) for i,index in enumerate(indices[0:len(indices)-1])] #returning annual maximum
    qmax.append(max(qq[indices[-1]:len(qq)]))
    return qmax


def get_best_fit_distribution(data,bins):
    """return a best fitted statistical distribution for peak flow; input: annual maximum Q (at least 10 valid points)"""
    if len(data)<10:
        print('... WARNING: sample size < 10, may not produce reliable distribution ...')
    
    # Get histogram of original data
    y, x = np.histogram(data, bins=bins, density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0
    
    # Best holders
    best_distribution = norm
    best_params = (0.0, 1.0)
    best_sse = np.inf  #sum of square error for choosing the best fit
    
    #finding the best distribution with the least sum of square error
    for distribution in distributions:
        # fit dist to data
        params = distribution.fit(data)
        
        # Separate parts of parameters
        arg = params[:-2]
        loc = params[-2]
        scale = params[-1]
        
        # Calculate fitted PDF and error with fit in distribution
        pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
        sse = np.sum(np.power(y - pdf, 2.0))
        
        # identify if this distribution is better
        if best_sse > sse > 0:
            best_distribution = distribution
            best_params = params
            best_sse = sse
    print(best_distribution.name)
    print(best_params)
    return (best_distribution.name, best_params)

def make_pdf(dist, params, size=5000):
    """Generate distributions's Probability Distribution Function """
    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)
    
    interval = float((end-start))/size
    return pdf,interval


def find_nearest_index(array, value):
    """Find nearest index in an array whose value is closest to input value"""
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx #return the index of the nearest valu


def get_bankfull_discharge(qq):
    """2-year return period flood to approximate bankfull Q; input data: daily Q, need to calculate fitted distribution's CDF"""
    qmax = get_annual_max_q(qq)
    bins = 50
    dist_name,params = get_best_fit_distribution(qmax,bins)
    best_dist = getattr(scipy.stats, dist_name)
    pdf, interval = make_pdf(best_dist,params)
    cdf = pdf.cumsum()*interval
    data = cdf.reset_index()
    data.columns = ['Q','CDF']
    data['PDF'] = pdf.values
    idx = find_nearest_index(data['CDF'].values,0.5)  #50% exceedance probab.; 2-yr return period flood to approximate bankfull
    q2 = data['Q'][idx]
    idx = find_nearest_index(data['CDF'].values,0.8)  #20% exceedance probab.; 5-yr return period flood to approximate bankfull
    q5 = data['Q'][idx]
    return (q2,q5)


if __name__ == '__main__':
    pfaf = 1
    print('...PFAF = %02d ...'%pfaf)
    fin = "chunked_pfaf_%02d.nc"%pfaf #netcdf chunking before (nccopy -c"time/100,COMID/1000" input.nc output.nc) for efficient slicing
    ds = xr.open_dataset(fin, chunks = {'time': 100, 'COMID':1000})
    comid = ds.COMID
    nf = len(comid)
    print('... nf = %s ...'%nf)

    Q2_list = []
    Q5_list = []
    QMEAN_list = []
    for i in range(nf)[rank::size]:
        comidnow = ds.COMID[i].values
        print('... COMID = %s ...'%comidnow)
#         import pdb;pdb.set_trace()
        qq = ds.Q[:,i].values
        q2,q5 = get_bankfull_discharge(qq)
        QMEAN_list.append(np.mean(qq))
        Q2_list.append(q2)
        Q5_list.append(q5)
#         df = pd.DataFrame({'QMEAN':[qmean],'Q2':[q2]})
#         df.to_csv('tmp/bankfull_q_%s.csv'%comidnow,index=False)

    df = pd.DataFrame({'COMID':comid[:nf].values,'QMEAN':QMEAN_list,'Q2':Q2_list,'Q5':Q5_list})
    df.to_csv('bankfull_q_pfaf_%02d.csv'%pfaf,index=False)
    
       
