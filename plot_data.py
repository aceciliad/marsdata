#!/usr/bin/env python
# -*- coding: utf-8 -*-

#==========================================================================
# Plot (overview) data
#==========================================================================

import numpy as np
import datetime, os, argparse, glob

from obspy import read
from obspy.signal.util import next_pow_2

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from matplotlib import cm
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from matplotlib import mlab as mlab

from get_times import *

def arguments():
    '''
    arguments
    '''
    ap = argparse.ArgumentParser(description='Plot event')

    ap.add_argument('--id', type=str, dest='id_data', help='id data', required=True)

    return ap.parse_args()


def detick(tr, detick_nfsamp,
           fill_val=None,
           freq_tick=1.0):
    # simplistic deticking by muting detick_nfsamp freqeuency samples around
    # 1Hz
    tr_out = tr.copy()
    Fs = tr.stats.sampling_rate
    NFFT = next_pow_2(tr.stats.npts)
    tr.detrend()
    df = np.fft.rfft(tr.data, n=NFFT)
    idx_1Hz = np.argmin(np.abs(np.fft.rfftfreq(NFFT) * Fs - freq_tick))
    if fill_val is None:
        fill_val = (df[idx_1Hz - detick_nfsamp - 1] + \
                    df[idx_1Hz + detick_nfsamp + 1]) / 2.
    df[idx_1Hz - detick_nfsamp:idx_1Hz + detick_nfsamp] /= \
        df[idx_1Hz - detick_nfsamp:idx_1Hz + detick_nfsamp] / fill_val
    tr_out.data = np.fft.irfft(df)[:tr.stats.npts]
    return tr_out


class COMPARE_WAVEFORMS():
    '''
    compare raw and deglitched waveforms
    '''

    def __init__(self, id_data):
        self.id_data= id_data
        self.get_raw()
        self.plot_data()
        return

    def get_raw(self,*,
                path_data='DATA'):

        '''
        get raw waveforms for chosen event
        '''
        try:
            # check if processed raw traces are available
            fname   = '*ZNE_{}.mseed'.format(self.id_data)
            files   = glob.glob(os.path.join(path_data,
                                             fname))
            if not files:
                raise ValueError('  Non processed deglitched streams found')
            if len(files)>1:
                raise ValueError(f'  More than 1 file for {self.id_data}')

            self.stream_raw = read(files[0])

        except:
            raise FileNotFoundError


        return


    def plot_data(self,*,
                  f0=1./500,
                  f1=10.,
                  dspace=1.5,
                  winlen=3600,
                  overlap=0.2):

        ax_st, ax_spe = self.set_figure(dspace=dspace)

        stream_raw = self.stream_raw.copy()
        stream_raw.filter('bandpass',
                          freqmin=f0, freqmax=f1,
                          corners=4, zerophase=True)

        kwargs_rw   = {'color':'dimgrey',
                       'lw':0.9,
                       'zorder':1}

        tref,_    = get_sol_start_end_utc(self.id_data)
        fact_time   = 1./(60*60)    # x axis in hours

        for ii, tr in enumerate(stream_raw):
            tr.normalize()
            ax_st.plot(tr.times(reftime=tref)*fact_time,
                       2*dspace-ii*dspace+tr.data, clip_on=True,
                       **kwargs_rw)
            ax_st.text(ax_st.get_xlim()[1], 2*dspace-ii*dspace+dspace/2,
                       tr.stats.component,
                       color='k', clip_on=False,
                       ha='right', va='top',
                       bbox= dict(boxstyle='square',
                                  facecolor='w', alpha=1,
                                  edgecolor='none'))


        stream_spec = self.stream_raw.copy()
        stream_spec.filter('bandpass',
                          freqmin=f0, freqmax=2*f1,
                          corners=4, zerophase=True)

        #stream_spec.resample(sampling_rate=2*f1, window='hann')

        tr = stream_spec.select(component='Z')[0]
        tr = detick(tr, detick_nfsamp=5)
        p, f, t = mlab.specgram(tr.data, NFFT=winlen,
                                Fs=tr.stats.sampling_rate,
                                noverlap=int(winlen*overlap),
                                pad_to=next_pow_2(winlen)*4)

        t -= (tref-tr.stats.starttime)

        flim    = np.array([0.001, 10])
        bol = np.array((f > flim[0], f < flim[1])).all(axis=0)

        ax_spe.pcolormesh(t/3600, f[bol], 10*np.log10(p[bol, :]),
                          norm=self.norm_spe,
                           cmap=self.cmap_spe,
                          rasterized=True)

        plt.savefig(f'Fig_sol{self.id_data}.pdf')
        return


    def set_figure(self,*,
                   fmin=0.01,
                   fmax=5,
                   dspace=1.5):

        fact_tr     = 2.e5
        fig = plt.figure()
        fig.set_size_inches(7,5)

        # Waveforms
        gs_st   = gridspec.GridSpec(ncols=1, nrows=3,
                                    bottom=0.1, top=0.96,
                                    left=0.09, right=0.84, figure=fig,
                                    wspace=0.1, hspace=.1)
        gs_cb   = gridspec.GridSpec(ncols=1, nrows=gs_st.nrows,
                                    bottom=gs_st.bottom,
                                    top=gs_st.top,
                                    left=gs_st.right+0.01,
                                    right=gs_st.right+0.03,
                                    figure=fig,
                                    hspace=gs_st.hspace)

        ax_st   = fig.add_subplot(gs_st[:2])
        ax_spe  = fig.add_subplot(gs_st[2])

        cb_spe  = fig.add_subplot(gs_cb[2])

        ax_spe.get_shared_x_axes().join(ax_spe, ax_st)

        ax_st.spines[['right', 'top', 'bottom']].set_visible(False)

        ax_st.set_ylim([-dspace/2, 2*dspace+dspace/2])
        ax_st.set_xlim([0, 24.7])
        ax_st.tick_params(labelbottom=False,
                          bottom=False,
                          labelleft=False,
                          left=False)
        ax_st.yaxis.set_major_locator(MultipleLocator(1.))
        #ax_st.yaxis.set_minor_locator(MultipleLocator(1000.))

        # Ax Spectrogram
        ax_spe.set_yscale('log')
        ax_spe.set_ylim([fmin, fmax])
        ax_spe.set_ylabel('Frequency [Hz]')
        ax_spe.set_xlabel('Time [hr]')
        ax_spe.xaxis.set_major_locator(MultipleLocator(6.))
        ax_spe.xaxis.set_minor_locator(MultipleLocator(1.))

        self.cmap_spe   = mpl.cm.plasma
        vmin_spe    = -210
        vmax_spe    = -120
        self.norm_spe    = mpl.colors.Normalize(vmin=vmin_spe,
                                                vmax=vmax_spe)
        cb = mpl.colorbar.ColorbarBase(cb_spe,
                                       orientation='vertical',
                                       cmap=self.cmap_spe,
                                       norm=self.norm_spe)
        cb.set_label(r'PSD [dB]')
        cb.ax.xaxis.set_major_locator(MultipleLocator(30))
        cb.ax.xaxis.set_minor_locator(MultipleLocator(10))

        return ax_st, ax_spe



if __name__=='__main__':

    results = arguments()
    id_data = results.id_data

    dir_figures = 'Figures'
    dir_data    = 'DATA'
    os.makedirs(dir_figures, exist_ok=True)

    # To download raw waveforms and get event info:
    COMPARE_WAVEFORMS(id_data)





