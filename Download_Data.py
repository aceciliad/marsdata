#!/usr/bin/env python
# -*- coding: utf-8 -*-

#==========================================================================
# Download data and generate input file to run SeisGlitch
#==========================================================================

import numpy as np
import datetime, os, argparse
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read, read_inventory

from get_times import *

import warnings
warnings.filterwarnings("ignore")

def arguments():
    '''
    arguments
    '''
    ap = argparse.ArgumentParser(description='Plot event')

    ap.add_argument('--id', type=str, dest='id_data', help='id data', required=True)

    return ap.parse_args()


def check_and_trim_traces(stream):
    # Check if the stream is empty
    if len(stream) == 0:
        return stream

    # Get the start and end time of the first trace
    reference_start_time = stream[0].stats.starttime
    reference_end_time = stream[0].stats.endtime

    # Iterate through the remaining traces in the stream
    for trace in stream[1:]:
        # Compare start time
        if trace.stats.starttime > reference_start_time:
            reference_start_time = trace.stats.starttime

        # Compare end time
        if trace.stats.endtime < reference_end_time:
            reference_end_time = trace.stats.endtime

    # Trim traces to the reference start and end time
    print(' Traces have different lenght.... Trimming')
    def trim_traces(stream, start_time, end_time):
        for trace in stream:
            trace.trim(starttime=start_time, endtime=end_time,
                       pad=True, fill_value=0.0)
            #trace.stats.endtime = trace.stats.starttime + trace.stats.npts / trace.stats.sampling_rate

        return stream

    stream  = trim_traces(stream, reference_start_time, reference_end_time)

    return stream


class DOWNLOAD_WAVEFORMS():
    '''
    download traces given a tstart and tend
    '''

    def __init__(self, id_data):
        self.id_data    = id_data
        self.download_data()
        return

    def read_timestamp(self,*,
                       hr_min=1,
                       hr_max=1,
                       data_path='DATA'):
        '''
        time stamps for downloading seismic data
        '''
        # start time to download
        # current time correspond to event S1222a and it will select data
        # from 0 hours before the event to 12 hours after its start
        #self.tref   = UTCDateTime(2022, 5, 4, 23, 27, 45, 836925)

        self.tref, self.tend = get_sol_start_end_utc(self.id_data)

        import ipdb; ipdb.set_trace()  # noqa
        self.starttime  = self.tref - datetime.timedelta(0, hours=hr_min)
        self.endtime    = self.tend + datetime.timedelta(0, hours=hr_max)

        return


    def download_data(self, *,
                      data_path='DATA',
                      type_st='VEL',
                      channel='BH*',
                      client_id='IRIS',
                      network='XB',
                      station='ELYSE',
                      location='02'):
        '''
        download data from IRIS and apply basic preprocessing
        '''
        self.read_timestamp()
        os.makedirs(data_path,
                    exist_ok=True)

        fnam_raw    = '{network}.{location}.{station}_UVW_'+self.id_data+'.mseed'
        fnam_pros   = '{network}.{location}.{station}_ZNE_'+self.id_data+'.mseed'

        print(f' >>> Downloading data {self.id_data}')
        # Set up client
        client  = Client(client_id)

        # retrieve stream
        self.stream  = client.get_waveforms(network=network,
                                            station=station,
                                            location=location,
                                            channel=channel,
                                            starttime=self.starttime,
                                            endtime=self.endtime,
                                            attach_response=True)
        if len(self.stream)>3:
            self.stream.merge()

        self.stream = check_and_trim_traces(self.stream)
        tstart  = self.stream[0].stats.starttime
        tend    = self.stream[0].stats.endtime

        for tr in self.stream:
            tr.stats.starttime  = tstart

        # Write raw traces to .mseed file
        fnam = os.path.join(data_path,
                            fnam_raw.format(**self.stream[0].stats))
        self.stream.write(fnam, format='mseed')

        # Traces are downloaded in the three original channels:
        # U, V, W (Longnonne et al., 2019)
        self.inventory = read_inventory(os.path.join(data_path,
                                                     'BH_inventory.xml'))

        self.stream.detrend(type='demean')
        self.stream.taper(0.1)
        self.stream.detrend()

        # Remove instrument response
        # use default values used by MQS
        print('     Removing instrument response')
        pre_filt    = [0.001, 0.005, 45, 50]
        self.stream.remove_response(pre_filt=pre_filt,
                                    output=type_st)
                                    #water_level=60)

        # Rotate to ZNE system
        print('     Rotating from UVW to ZNE')
        self.stream.rotate('->ZNE',
                           components=['UVW'],
                           inventory=self.inventory)

        # Write pre-processed traces to .mseed file
        fnam = os.path.join(data_path,
                            fnam_pros.format(**self.stream[0].stats))
        self.stream.write(fnam, format='mseed')

        return


if __name__=='__main__':

    results = arguments()
    id_data = results.id_data

    #list_sols   = np.arange(184, 200, 1)
    #sols    = [str(xx) for xx in list_sols]

    # To download raw waveforms and get event info:
    DOWNLOAD_WAVEFORMS(id_data)



