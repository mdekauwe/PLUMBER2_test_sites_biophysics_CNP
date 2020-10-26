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
import datetime
import os
import glob
from optparse import OptionParser

def main(site, met_fname, flux_fname, biophysics_fname, C_fname, CN_fname,
         CNP_fname):

    df1 = read_cable_file(biophysics_fname, type="CABLE")
    df1 = resample_timestep(df1, type="CABLE")

    df2 = read_cable_file(C_fname, type="CABLE")
    df2 = resample_timestep(df2, type="CABLE")

    df3 = read_cable_file(CN_fname, type="CABLE")
    df3 = resample_timestep(df3, type="CABLE")

    df4 = read_cable_file(CNP_fname, type="CABLE")
    df4 = resample_timestep(df4, type="CABLE")

    df_flx = read_cable_file(flux_fname, type="FLUX")
    df_flx = resample_timestep(df_flx, type="FLUX")

    #df_met = read_cable_file(met_fname, type="MET")
    #df_met = resample_timestep(df_met, type="MET")

    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(hspace=0.05)
    fig.subplots_adjust(wspace=0.2)
    plt.rcParams['text.usetex'] = False
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.sans-serif'] = "Helvetica"
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['font.size'] = 12
    plt.rcParams['legend.fontsize'] = 12
    plt.rcParams['xtick.labelsize'] = 12
    plt.rcParams['ytick.labelsize'] = 12


    colours = plt.cm.Set2(np.linspace(0, 1, 7))

    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    #axx1 = ax1.twinx()
    #axx2 = ax2.twinx()


    axes = [ax1, ax2,]
    #axes2 = [axx1, axx2,]
    vars = ["GPP", "Qle"]
    #for a, x, v in zip(axes, axes2, vars):
    for a, v in zip(axes, vars):
        if v == "Qle":
            a.plot(df_flx[v].index.to_pydatetime(),
                   df_flx[v].rolling(window=7).mean(), c="seagreen", lw=2.0,
                   ls="-", label="Observations")
        a.plot(df1[v].index.to_pydatetime(), df1[v], c="orange",
               lw=1.5, ls="-", label="Biophysics")
        a.plot(df2[v].index.to_pydatetime(), df2[v]), c="lightblue",
               lw=1.5, ls="-", label="C")
        a.plot(df3[v].index.to_pydatetime(), df3[v], c="DodgerBlue",
               lw=1.5, ls="-", label="CN")
        a.plot(df4[v].index.to_pydatetime(), df4[v], c="darkblue",
               lw=1.5, ls="-", label="CNP")

        #x.bar(df_met.index, df_met["Rainf"], alpha=0.3, color="black")
    #ax2.set_ylim(0, 200)
    labels = ["GPP (g C m$^{-2}$ d$^{-1}$)", "LE (W m$^{-2}$)"]
    for a, l in zip(axes, labels):
        a.set_ylabel(l, fontsize=12)

    #for x, l in zip(axes2, labels):
    #    x.set_ylabel("Rainfall (mm d$^{-1}$)", fontsize=12)


    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2.legend(numpoints=1, loc="best", ncol=2, fontsize=10)


    for a in axes:
        #a.set_xlim([datetime.date(2002,10,1), datetime.date(2003, 4, 1)])
        #a.set_xlim([datetime.date(2002,12,1), datetime.date(2003, 5, 1)])

        a.set_xlim([datetime.date(2008,1,1), datetime.date(2009, 1, 1)])
        #a.set_xlim([datetime.date(2006,11,1), datetime.date(2007, 4, 1)])

    plot_fname = "/Users/mdekauwe/Desktop/CNP_timeseries.png"
    if plot_fname is None:
        plt.show()
    else:
        #fig.autofmt_xdate()
        fig.savefig(plot_fname, dpi=150, bbox_inches='tight', pad_inches=0.1)


def read_cable_file(fname, type=None):

    if type == "CABLE":
        vars_to_keep = ['GPP','Qle','LAI','TVeg', 'ESoil','NEE']
    elif type == "FLUX":
        vars_to_keep = ['GPP','Qle']
    elif type == "MET":
        vars_to_keep = ['Rainf']

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

    if type == "CABLE":
        ds['TVeg'] *= float(time_jump)
        ds['ESoil'] *= float(time_jump)
        ds['GPP'] *= float(time_jump)
    elif type == "FLUX":
        ds['GPP'] *= float(time_jump)
    elif type == "MET":
        ds['Rainf'] *= float(time_jump)


    df = ds.to_dataframe()

    start = reference_date.strip().split(" ")[0].replace("-","/")
    df['dates'] = pd.date_range(start=start, periods=len(df), freq=freq)
    df = df.set_index('dates')

    return df


def resample_timestep(df, type=None):

    UMOL_TO_MOL = 1E-6
    MOL_C_TO_GRAMS_C = 12.0


    if type == "CABLE":
        # umol/m2/timestep -> g/C/timestep
        df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C


        method = {'GPP':'sum', 'TVeg':'sum', "Qle":"mean"}
    elif type == "FLUX":
        # umol/m2/timestep -> g/C/timestep
        df['GPP'] *= UMOL_TO_MOL * MOL_C_TO_GRAMS_C
        method = {'GPP':'sum', "Qle":"mean"}


    df = df.resample("D").agg(method)

    return df


    return dates

if __name__ == "__main__":

    #site = "TumbaFluxnet"
    output_dir = "outputs"
    #flux_dir = "../../flux_files/"
    #met_dir = "../../met_data/plumber_met/"
    met_dir = "/Users/mdekauwe/research/OzFlux"
    flux_dir = "/Users/mdekauwe/research/OzFlux"

    """
    parser = OptionParser()
    parser.add_option("-a", "--fname1", dest="fname1",
                      action="store", help="filename",
                      type="string",
                      default=os.path.join(output_dir, "original_out.nc"))
    parser.add_option("-b", "--fname2", dest="fname2",
                      action="store", help="filename",
                      type="string",
                      default=os.path.join(output_dir, "%s_out.nc" % site))
    parser.add_option("-c", "--fname3", dest="fname3",
                      action="store", help="filename",
                      type="string",
                      default=os.path.join(output_dir, "%s_out.nc" % site))
    parser.add_option("-p", "--plot_fname", dest="plot_fname", action="store",
                      help="Benchmark plot filename", type="string")
    (options, args) = parser.parse_args()
    """

    flux_fname = os.path.join(flux_dir, "TumbarumbaOzFlux2.0_flux.nc")
    met_fname = os.path.join(met_dir, "TumbarumbaOzFlux2.0_met.nc")


    site = "AU-Tum"
    biophysics_fname = "outputs/AU-Tum_2002-2016_OzFlux_Met_out.nc"
    C_fname = "outputs/%s_C_out_cable_simulation.nc" % (site)
    C_fname = "outputs/%s_C_out_cable_simulation.nc" % (site)
    CN_fname = "outputs/%s_CN_out_cable_simulation.nc" % (site)
    CNP_fname = "outputs/%s_CNP_out_cable_simulation.nc" % (site)

    main(site, met_fname, flux_fname, biophysics_fname, C_fname, CN_fname,
         CNP_fname)
