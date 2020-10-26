#!/usr/bin/env python

"""
Get casa values

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (21.07.2018)"
__email__ = "mdekauwe@gmail.com"

import matplotlib.pyplot as plt
import sys
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator
import os
import xarray as xr
import glob

def open_casa_and_add_time(fname, start_date):
    ds = xr.open_dataset(fname)
    N = len(ds.Nupland)
    first_date = pd.to_datetime(start_date)
    time = first_date + pd.to_timedelta(np.arange(N), 'D')
    ds['time'] = time

    return (ds)

if __name__ == "__main__":

    cycle = "CN"
    fname = "*_%s_out_casa_simulation.nc" % (cycle)
    fname = glob.glob(os.path.join("outputs", fname))[0]
    ds = open_casa_and_add_time(fname, start_date="01/01/2002")
    #fname = "*_%s_out_casa_historical.nc" % (cycle)
    #fname = glob.glob(os.path.join("outputs", fname))[0]
    #ds = open_casa_and_add_time(fname, start_date="01/01/1850")


    cf = ds.cplant[:,0].groupby("time.year").mean().values
    cw = ds.cplant[:,1].groupby("time.year").mean().values
    cr = ds.cplant[:,2].groupby("time.year").mean().values
    ncf = ds.nplant[:,0].groupby("time.year").mean().values
    ncw = ds.nplant[:,1].groupby("time.year").mean().values
    ncr = ds.nplant[:,2].groupby("time.year").mean().values

    print("\n\nPlant C & N")
    print(np.mean(cf), np.mean(cw), np.mean(cr), np.mean(cw/cf))
    print(np.mean(ncf), np.mean(ncw), np.mean(ncr))
    print(np.mean(ncf)/np.mean(cf), np.mean(ncw)/np.mean(cw),
          np.mean(ncr)/np.mean(cr))


    cactive = ds.csoil[:,0].groupby("time.year").mean().values
    cslow = ds.csoil[:,1].groupby("time.year").mean().values
    cpassive = ds.csoil[:,2].groupby("time.year").mean().values
    nactive = ds.nsoil[:,0].groupby("time.year").mean().values
    nslow = ds.nsoil[:,1].groupby("time.year").mean().values
    npassive = ds.nsoil[:,2].groupby("time.year").mean().values

    print("\nSoil C & N")
    print(np.mean(cactive), np.mean(cslow), np.mean(cpassive))
    print(np.mean(nactive), np.mean(nslow), np.mean(npassive))
    print(np.mean(nactive)/np.mean(cactive), np.mean(nslow)/np.mean(cslow),
          np.mean(npassive)/np.mean(cpassive))


    nmin = ds.Nsnet[:].groupby("time.year").sum().values
    ndep = ds.Nmindep[:].groupby("time.year").sum().values
    nfix = ds.Nminfix[:].groupby("time.year").sum().values
    input = nmin + nfix + ndep

    nleach = ds.Nminleach[:].groupby("time.year").sum().values
    nup = ds.Nupland[:].groupby("time.year").sum().values
    ndep = ds.Nmindep[:].groupby("time.year").sum().values
    nloss = ds.Nminloss[:].groupby("time.year").sum().values
    loss = nloss + nleach + nup

    print("\nChange in NpoolM")
    print(np.mean(input), np.mean(loss), np.mean(input-loss))

    nimmobilisation = ds.Nsimm[:].groupby("time.year").mean().values
    print("\nN immobilisation")
    print(np.mean(nimmobilisation))
    print("\nN Uptake - GDAY Eucf 10 g N m-2 y-1")
    print(np.mean(nup))
    print("\nN gross mineralisation - Eucf 10 g N m-2 y-1")
    print(np.mean(nup))
