#!/usr/bin/env python

"""
Plot visual benchmark (average seasonal cycle) of old vs new model runs.

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (18.10.2017)"
__email__ = "mdekauwe@gmail.com"

import xarray as xr
import matplotlib.pyplot as plt
import sys
import datetime as dt
import pandas as pd
import numpy as np
from matplotlib.ticker import FixedLocator
import os

def main(biophysics_fname, C_fname, CN_fname=None, CNP_fname=None,
         plot_fname=None):

    df_biophys = read_cable_file(biophysics_fname)
    df_biophys = resample_to_seasonal_cycle(df_biophys)

    df_C = read_cable_file(C_fname)
    df_C = resample_to_seasonal_cycle(df_C)
    if CN_fname is not None:
        df_CN = read_cable_file(CN_fname)
        df_CN = resample_to_seasonal_cycle(df_CN)
    if CNP_fname is not None:
        df_CNP = read_cable_file(CNP_fname)
        df_CNP = resample_to_seasonal_cycle(df_CNP)

    fig = plt.figure(figsize=(10,9))
    fig.subplots_adjust(hspace=0.3)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12

    ax1 = fig.add_subplot(3,2,1)
    ax2 = fig.add_subplot(3,2,2)
    ax3 = fig.add_subplot(3,2,3)
    ax4 = fig.add_subplot(3,2,4)
    ax5 = fig.add_subplot(3,2,5)
    ax6 = fig.add_subplot(3,2,6)

    axes = [ax1, ax2, ax3, ax4, ax5, ax6]
    vars = ["GPP", "NEE", "Qle", "LAI", "TVeg", "ESoil"]
    for a, v in zip(axes, vars):
        a.plot(df_biophys.month, df_biophys[v], c="orange", lw=2.0, ls="-",
               label="Biophysics")
        a.plot(df_C.month, df_C[v], c="lightblue", lw=2.0, ls="-", label="C")
        a.plot(df_CN.month, df_CN[v], c="DodgerBlue", lw=2.0, ls="-",
               label="CN")
        a.plot(df_CNP.month, df_CNP[v], c="darkblue", lw=2.0, ls="-",
               label="CNP")

    labels = ["GPP (g C m$^{-2}$ d$^{-1}$)", "NEE (g C m$^{-2}$ d$^{-1}$)",\
              "Qle (W m$^{-2}$)", "LAI (m$^{2}$ m$^{-2}$)",\
              "TVeg (mm d$^{-1}$)", "Esoil (mm d$^{-1}$)"]
    for a, l in zip(axes, labels):
        a.set_title(l, fontsize=12)

    xtickagaes_minor = FixedLocator([2, 3, 4, 5, 7, 8, 9, 10, 11])
    for i,a in enumerate(axes):
        a.set_xticks([1, 6, 12])
        if i != 1 and i != 2 and i != 4 and i != 5:
            a.set_ylim(ymin=0)
        a.xaxis.set_minor_locator(xtickagaes_minor)
        a.set_xticklabels(['Jan', 'Jun', 'Dec'])
        if i < 4:
            plt.setp(a.get_xticklabels(), visible=False)
    ax2.legend(numpoints=1, loc="best", fontsize=10, ncol=2)



    plot_dir = "plots"
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    fig.savefig(os.path.join(plot_dir, plot_fname), bbox_inches='tight',
                pad_inches=0.1)

def read_cable_file(fname):

    vars_to_keep = ['GPP','Qle','LAI','TVeg', 'ESoil','NEE']

    ds = xr.open_dataset(fname, decode_times=False)

    time_jump = int(ds.time[1].values) - int(ds.time[0].values)

    if time_jump == 3600:
        freq = "H"
    elif time_jump == 1800:
        freq = "30M"
    elif time_jump == 86400:
        freq = "D"
    else:
        raise("Time problem")

    units, reference_date = ds.time.attrs['units'].split('since')
    ds = ds[vars_to_keep].squeeze(dim=["x","y"], drop=True)

    #ds['TVeg'] *= float(time_jump)
    #ds['ESoil'] *= float(time_jump)
    #ds['GPP'] *= float(time_jump)

    df = ds.to_dataframe()

    start = reference_date.strip().split(" ")[0].replace("-","/")
    df['dates'] = pd.date_range(start=start, periods=len(df), freq=freq)
    df = df.set_index('dates')

    return df



def resample_to_seasonal_cycle(df, OBS=False):

    UMOL_TO_MOL = 1E-6
    MOL_C_TO_GRAMS_C = 12.0
    SEC_2_DAY = 86400.

    # umol/m2/s -> g/C/d
    df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_DAY
    df['NEE'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_DAY

    # kg/m2/s -> mm/d
    df['TVeg'] *= SEC_2_DAY
    df['ESoil'] *= SEC_2_DAY

    method = {'GPP':'mean', 'NEE':'mean', 'Qle':'mean', 'LAI':'mean',
              'TVeg':'mean', 'ESoil':'mean'}
    df = df.resample("M").agg(method).groupby(lambda x: x.month).mean()
    df['month'] = np.arange(1,13)

    return df

if __name__ == "__main__":

    site = "AU-Tum"
    biophysics_fname = "outputs/%s_biophysics.nc" % (site)
    C_fname = "outputs/%s_C_out_cable_simulation.nc" % (site)
    C_fname = "outputs/%s_C_out_cable_simulation.nc" % (site)
    CN_fname = "outputs/%s_CN_out_cable_simulation.nc" % (site)
    CNP_fname = "outputs/%s_CNP_out_cable_simulation.nc" % (site)
    plot_fname = "%s_seasonal_plot_C_vs_CN_vs_CNP.png"  % (site)
    main(biophysics_fname, C_fname, CN_fname=CN_fname, CNP_fname=CNP_fname,
         plot_fname=plot_fname)
