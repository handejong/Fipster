#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Welcome to Fipster Python. The most important part of this FIP analysis code is the
FIP_signal python class. You instantiate one like so:

>>> signal = FIP_signal(filename = 'your_filename.mat')

To get an overview of the available methods type:

>>> help(FIP_signal)

Created on Tue Nov 12 14:41:27 2019

@author: Han de Jong (j.w.dejong@berkeley.edu)
"""

# These are neccesary for file browsing
from scipy.io import loadmat
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import os
import pandas as pd
import seaborn as sns
from scipy import signal as scipy_signal
from sklearn import metrics
import warnings


class FIP_signal:
    # The FIP_signal object contains 1 or multiple FIP signals.
    
    def __init__(self, **kwargs):
        # Constructor
        
        # Set basic properties
        self.hasdata = False
        self.filepath = []
        self.nr_signals = 0
        self.signal = []
        self.raw_ref =[]
        self.framerate = []
        self.labels = []
        self.filename = []
        self.notes = []
        self.logAI = []
        self.detrend = False
        self.external_signal = []
        self.smooth = []

        # A dictionary that will contain all the timestamps
        self.timestamps = {}

        # This is where you can put a filter that would be applied
        self.filter = None

        # This is just for formatting of the figure
        self.facecolor = 'w'
        self.axcolor = 'k'
        self.colormap = colormaps['tab20'](np.linspace(0, 1, 20))
        
        # Set the settings. (IMPORTANT)
        self.settings = {
                'fit ref': 'polyfit',
                'rolling window': 1000,
                'signal unit': '∆F/F',}
   
        # Deal with input arguments
        for key, value in kwargs.items():
            
            if key=='filename':
                # open that file
                self.filepath = value
            elif key=='detrend':
                self.detrend = value
            else:
                print('unknown argument: {0}, ignored'.format(key))
                
        # Do we have an idea about what to import? -> Load that data
        if self.filepath != []:
            self.load_data(self.filepath)
                    
        # Is data imported? If not, open a file browser
        if not self.hasdata:
            print('No data loaded, created empty FIP_signal object.')
            #TODO include filebrowser
            
                    
    def load_data(self, filename):
        # Load .mat file into the object
        
        print('Loading: {0}'.format(filename))
        
        # Try two different ways of importing old and new Matlab files
        with h5py.File(filename, mode='r') as data:
              
            # Unpack individual variables  
            for key, value in data.items():
                
                if key == 'framerate':
                    self.framerate = int(value[()])
                
                if key == 'notes':
                    try:
                        string_value = [i[0] for i in data[value[0][0]][()]]
                        self.notes = ''.join(map(lambda x: chr(x), string_value))
                    except:
                        self.notes = ''
                
                if key == 'signal':
                    self.signal = np.array(value[()])
                
                if key == 'ref':
                    self.raw_ref = np.array(value[()])
                    
                if key == 'labels':
                    self.labels = []
                    for i in range(value.shape[1]):
                        string_value = [j[0] for j in data[value[0][i]][()]]
                        label = ''.join(map(lambda x: chr(x), string_value))
                        self.labels.append(label)
            
        # If there is a ref, divide the framerate by 2
        if not len(self.raw_ref)==0:
            self.framerate = self.framerate/2
        else: # If there is no ref, put some zeros and prevent ref fit
            self.raw_ref = np.zeros(self.signal.shape)
            self.settings['fit ref'] = 'no fit'
        
        # Make sure self.signal and self.ref have 3dims
        if len(self.signal.shape) == 2:
            self.signal = np.expand_dims(self.signal, axis = 2)
            self.raw_ref = np.expand_dims(self.raw_ref, axis = 2)
        else:
            self.signal = np.swapaxes(self.signal, 0, 1)
            self.signal = np.swapaxes(self.signal, 1, 2)
            self.raw_ref = np.swapaxes(self.raw_ref, 0, 1)
            self.raw_ref = np.swapaxes(self.raw_ref, 1, 2)
            
        # Is there an associated logAI file?
        new_filename = filename[0:-4] + '_logAI.csv'
        if os.path.isfile(new_filename):
            self.load_logAI(new_filename)
            
        # Finish up
        self.nr_signals = self.signal.shape[2]
        self.hasdata = True
        self.filename = filename
        self.external_signal = [False] * self.signal.shape[2]
        self.smooth = [False] * self.signal.shape[2]

        # Print recording notes
        print(' ')
        print(f'RECORDING NOTES: {self.notes}')
        print(' ')


    def load_logAI(self, filename):
        """
        Load the associated logAI file with the FIP signal.

        The logAI file contains analog input signals collected on the NI Daq during
        a recording. It is generally sampled at 200Hz (but this can vary). There are
        some quircks in how Matlab stores timepoints at which logAI data is collected
        which is annoying and the reason we have to resample.

        """
        print(f'Loading: {filename}')
        logAI = pd.read_csv(filename)
        
        # Name the columns
        logAI.columns = ['Time'] + [f'TTL {i}' for i in range(1, logAI.shape[1])] 
        
        # Remove columns if there is no signal above 1mV
        indexer = logAI.max(axis=0)>1
        logAI = logAI.loc[:, indexer].copy()

        # How about that timeline, can we safely resample?
        sampling_interval = round(logAI['Time'].iloc[-1]/len(logAI['Time']), 7)
        print('Sampling interval is: {0}ms'.format(sampling_interval*1000))
        if sampling_interval == 0.005000:
            logAI.Time = np.linspace(0.005, len(logAI)*0.005, len(logAI))
            print('...resampled')
        elif (sampling_interval>0.0049) & (sampling_interval<0.0051):
            logAI.time = np.linspace(0.005, len(logAI)*0.005, len(logAI))
            print('...resampled ANYWAY')

        self.logAI = logAI


    def import_signal(self, new_signal):
        # Import external signal
        
        # Error handeling
        # TODO
        #   .... does it have it's own timeline?
        #   .... is the sample rate equal to the internal signal
        #   .... does it come with a reference?
        if not len(new_signal) == self.signal.shape[1]:
            print("This function is not finished")
            return 0
        
        # Include the new signal
        temp = np.zeros([2, self.signal.shape[1], self.nr_signals+1])
        temp[:, :, :self.nr_signals] = self.signal
        temp[1, :, self.nr_signals] = temp[1, :, 0]
        temp[0, :, self.nr_signals] = new_signal
        self.signal = temp
        
        # Include the new ref (currently just leave it 0)
        temp = np.zeros([2, self.raw_ref.shape[1], self.nr_signals+1])
        temp[:, :, :self.nr_signals] = self.raw_ref
        temp[1, :, self.nr_signals] = temp[1, :, 0]
        self.raw_ref = temp
        
        # Finish up
        self.nr_signals += 1
        self.smooth.append(0)
        self.external_signal.append(True)
        # NOTE: set external_signal to True if you want to normalized it in the
        # same way as the internal signals.
        
         
    def check_settings(self):
        # Checks if all settings are valid (to prevent errors).
        
        # Cycle through the keys and check them all
        for key, value in self.settings.items():
            
            found_key = False
            
            if key == 'fit ref':
                found_key = True;
                options = ['polyfit', 'means', 'no fit', 'rolling polyfit']
                if not (value in options):
                    raise Exception("{0} is not a valid fit method. pick from: " \
                                    "\'polyfit\', \'means\', \'rolling polyfit\', or \'no fit\'.".format(value))
                    
            if key == 'signal unit':
                found_key = True;
                if not (value == '∆F/F' or value == 'Z-score' or value == 'deltaF'):
                    raise Exception("{0} is not a valid fit method. pick from: " \
                                    "\'∆F/F\', \'deltaF\', or \'Z-score\'.".format(value))
             
            if key == 'rolling window':
                found_key = True
                if not (value>1):
                    print('WARNING: a rolling of <1 is interpreted as a fraction of total.')
                        
            if not found_key:
                self.settings.pop(key)
                raise Exception("{0} is not a valid setting... removed".format(key))
  
    
    def get_data(self):
        # Calculate the data based on the signal, the ref and the settings
        
        # Check settings
        self.check_settings()
        
        # Calculate the normalized data based on the settings
        data = self.signal.copy()
        ref = self.get_ref()

        # Do we filter?
        if not self.filter is None:
            print('Applying provided filter to signal.')
            data[0, :, :] = scipy_signal.filtfilt(self.filter[0], self.filter[1], data[0, :, :], axis=0)
        
        # Subtract the fitted reference, add the mean value
        for i in range(0, self.nr_signals):

            if not self.settings['fit ref'] == 'no fit':
                data[0, :, i] = data[0, :, i] - ref[0, :, i] + np.mean(ref[0, :, i])
            
            # should we detrend? (Apply Savitzky-Golay filter)
            if self.detrend:
                data[0, :, i] = data[0, :, i] - \
                scipy_signal.savgol_filter(data[0, :, i], 5001, 1) + \
                np.min(data[0, :, i])
                
            # Do we smooth the signal
            if self.smooth[i]:
                data[0, :, i] = scipy_signal.savgol_filter(data[0, :, i], 
                                                     self.smooth[i], 1)
                
            # Do we want Z-scores?
            if self.settings['signal unit'] == 'Z-score':
                data[0, :, i] = (data[0, :, i] - np.mean(data[0, :, i]))/\
                np.std(data[0, :, i])
                
            # Report delta/F over F
            # We are assuming that we don't want to do this for external signals
            if not self.external_signal[i]:
                F = np.median(data[0, :, i])
                data[0, :, i] = (data[0, :, i]-F)/F
            
        # Return the normalized data
        return data
    
    
    def get_ref(self, detrend=False):
        """
        Calculate the fitted reference signal
        Normally you'd never want do de-trend (e.g. remove bleaching) from the
        reference trace, however it might be that you want to make peri-
        event plots of the reference trace in which case you'll want to
        remove the bleaching trend using detrend=True (the peri-event method
        does this automatically).
        """
        
        # Check settings
        self.check_settings()
        
        # Grab raw data
        ref = self.raw_ref.copy()
        signal = self.signal.copy()

        # Do we filter?
        if not self.filter is None:
            ref[0, :, :] = scipy_signal.filtfilt(self.filter[0], self.filter[1], ref[0, :, :], axis=0)
            signal[0, :, :] = scipy_signal.filtfilt(self.filter[0], self.filter[1], signal[0, :, :], axis=0)
        
        # Do we de-trend?(Apply Savitzky-Golay filter)        
        if detrend:
            window = min((5000, ref.shape[1]))
            for i in range(0, self.nr_signals):
                ref[0, :, i] = ref[0, :, i] - \
                scipy_signal.savgol_filter(ref[0, :, i], window, 1) + \
                np.min(ref[0, :, i])
               
        # Go over the 4 methods
        if self.settings['fit ref'] == 'polyfit':
            for i in range(0, self.nr_signals):
                p = np.polyfit(ref[0,:,i], signal[0, :, i], 1)
                ref[0, :, i] = ref[0, :, i] * p[0] + p[1]
                
        if self.settings['fit ref'] == 'means':
            for i in range(0, self.nr_signals):
                correction = np.mean(signal[0, :, i])/np.mean(ref[0, :, i])
                ref[0, :, i] = ref[0, :, i] = ref[0, :, i] * correction
        
        if self.settings['fit ref'] =='no fit':
            ref = ref
           
        if self.settings['fit ref'] == 'rolling polyfit':
            ref = self.rolling_polyfit()
        
        return ref   
              
    def rolling_polyfit(self):
        """
        This is just a 1st order polynomal fit, except we execute it on a
        rolling interval. This is good for longer recordings or when there
        is a shift in signal intensity between the signal and the refference
        """
        
        ref = self.raw_ref.copy()
        signal = self.signal

        # Make the output
        output = ref.copy()
        
        # Apply the fit over an interval of length 'rolling window'
        window = int(self.settings['rolling window']/2)
        for i in range(0, self.nr_signals): # For every signal
            for j in range(ref.shape[1]):
                start = max([j-window, 0])
                end = min([j+window, ref.shape[1]-1])
                p = np.polyfit(ref[0,start:end,i], 
                               signal[0, start:end, i], 1)
                output[0, j, i] = ref[0, j, i] * p[0] + p[1]
        
        return output   


    def plot(self, signals:list = 'all', raw_data = True, timestamps = True):
        """
        This is the main function to quickly plot all the available data in a
        Matplotlib figure.

        Parameters:
        -----------
        signals : list, tuple or string (if 'all')
            A list of all signals you want to plot.
        raw_data: Bool
            If True, this will plot raw non-normalized data. It will fit the
            reference to the signal though. You set that separately in "settings".
        timestamps: Bool
            You can plot either raw log_AI traces or timestamps. If set to True
            this will plot timestamps in the bottom plot.
        """
        
        # Error handling on the input arguments
        if len(self.timestamps)==0:
            timestamps = False
        
        # Figure out how many signals to show
        if signals=='all':
            signals = range(0, self.nr_signals)

        # Make a figure
        figure, axs = plt.subplots(len(signals)+1, 1, figsize = (16, 3*len(signals)+3), 
            tight_layout=True, sharex=True, facecolor=self.facecolor)
        _ = [self._default_formatting(ax) for ax in axs]
        
        # Plot the signals
        if raw_data:
            ref = self.get_ref()

            # Grab the signals
            to_plot = self.signal.copy()

            # Filter if requested
            if not self.filter is None:
                print('Applying provided filter to signal')
                to_plot[0, :, :] = scipy_signal.filtfilt(self.filter[0], self.filter[1], to_plot[0, :, :], axis=0)

            for i, signal_i in enumerate(signals):

                axs[i].plot(to_plot[1, :, signal_i], to_plot[0, :, signal_i], lw = 0.5,
                    label = 'signal')
                if not ref[0, :, signal_i].sum() == 0:
                    axs[i].plot(ref[1, :, signal_i], ref[0, :, signal_i], lw = 0.5, 
                        label='ref')
                axs[i].set_ylabel(f'{self.labels[signal_i]} (F)')
        
        # Plot normalized data instead
        else:
            data = self.get_data()
            for i, signal_i in enumerate(signals):
                axs[i].plot(data[1, :, signal_i], data[0, :, signal_i], lw = 0.5, label=f'Channel_{i}')
                axs[i].set_ylabel(f"{self.labels[signal_i]} ({self.settings['signal unit']})")

        # Some formatting
        axs[i].legend(loc = 'upper right')

        # Plot the logAI data
        if not timestamps:
            for i in range(1, self.logAI.shape[1]):
                axs[-1].plot(self.logAI.iloc[:, 0], self.logAI.iloc[:, i],\
                         label = self.logAI.columns[i], lw=0.5)
            axs[-1].set_ylabel('AI (mV)')    
        else:
            offset = 0
            for key, value in self.timestamps.items():
                axs[-1].eventplot(value, lineoffsets=offset, label=key,
                    colors=self.colormap[offset, :3])
                offset += 1
            axs[-1].set_ylabel('Stamp #')    
        axs[-1].set_xlabel('Time (s)')
        axs[-1].legend(loc = 'upper right')

        # Some formatting
        figure.suptitle(self.filename)

        return figure, axs

    
    def derive_timestamps(self, TTL = 'all', names = 'default', min_length = 0,
                          max_length = 9999, threshold = 1, store_stamps:bool = True):
        """
        DERIVE_TIMESTAMPS converts log AI voltage traces to discrete timestamps.

        Parameters:
        -----------
        TTL: Str, int or list of either
            A list or individual TTL channel.
        names: Str or list of Str
            A list of names to be used for the imported pulses
        min_length: numeric
            Minimum pulse length (s) to be included. Shorter pulses are ignored.
        max_length: numeric
            Maximum pulse length (s) to be included. Longer pulses are ignored.
        threshold: numeric
            Threshold (voltage) above which a pulse onset is detected.
        store_stamps: Bool
            Controls if the stamps are saved into the FIP_SIGNAL object.

        Returns:
        --------
        A dictionary with the names as keys and DataFrames containing the timestamps
        as values.
        """
        
        # Output variable
        output_stamps = {}

        # Individual columns or all columns
        if TTL == 'all':
            columns = self.logAI.columns[1:]
        elif TTL.__class__ == str:
            columns = [TTL]
        elif TTL.__class__ == int:
            columns = [f'TTL {TTL}']
        else:
            if TTL[0].__class__ == int:
                columns = [f'TTL {i}' for i in TTL]
            else:
                columns = TTL    
            
        # Figure out TTL stamps names
        if names=='default':
            names = columns
            
        # Just like with columns, make sure it's a list
        if names.__class__ == str:
            names = [names]
            
        # Error handling
        if not len(columns)==len(names):
            raise Exception('Number of timestamp series and names are unequal.')
            
        # Go through all the log AI
        timeline = self.logAI['Time'].values
        for i in range(0, len(columns)):
            stamps = {'start': [], 'end': [], 'length': []}
            
            col = columns[i]
            name = names[i]
            signal = self.logAI[col].values
            
            # Grab TTL_high timepoins
            TTL_high = signal>threshold
            dif = np.append([0], np.diff(TTL_high.astype(int)))
            
            # Grab the starts and ends
            starts = timeline[dif==1]
            ends = timeline[dif==-1]

            # Are there any stamps in here?
            if (len(starts)==0) | (len(ends)==0):
                print(f'No stamps on {columns[i]} ({names[i]})')
                continue
            
            # Error handling
            if ends[0]<starts[0]:
                ends = ends[1:]
            if starts[-1] > ends[-1]:
                starts = starts[:-1]
            if not len(starts) == len(ends):
                warnings.warn(f'For stamps {name} there is an unequal number of TTL starts and stops, this is an error!')
                print('Trying to match stamp starts with stamp ends')
                new_starts = []
                new_ends = []
                for stamp in starts:
                    try:
                        new_end = ends[ends>=stamp][0]
                        new_starts.append(stamp)
                        new_ends.append(new_end)
                    except:
                        pass
                starts = np.array(new_starts)
                ends = np.array(new_ends)
            
            # Make the set of stamps
            stamps['start'] = starts
            stamps['end'] = ends
            stamps['TTL'] = col
            stamps['name'] = name
            if store_stamps:
                self.timestamps[name+'_onset'] = starts
                self.timestamps[name+'_offset'] = ends
                print('Importing {0} stamps, named: {1}'.format(len(stamps['start']), stamps['name']))
            else:
                print(f"Derived {len(stamps['start'])} stamps with name {stamps['name']}.")
            
            # Can we calculate the stamp length?
            stamps['length'] = stamps['end'] - stamps['start']
            selector = stamps['length']>=min_length
            selector[stamps['length']>max_length] = False
            
            # Apply selector
            stamps['start'] = stamps['start'][selector]
            stamps['end'] = stamps['end'][selector]
            stamps['length'] = stamps['length'][selector]

            # Add these stamps to the output
            output_stamps[name] = pd.DataFrame(stamps)
            
        return output_stamps


    def quick_peri(self, TTL='all', window = 10, from_ref = False,
        stamps_min = 0, stamps_max = np.inf, offset = False):
        """
        QUICK_PERI will make a quick peri-event plot around the pulses
        captures in TTL. Here TLL can be an int or a string.

        Example:

            >>> PE = signal.quick_peri(1)

        This will make the sweepset object around the TTL stamps on
        TTL 1. It will them plot this dataset.

        Parameters:
        -----------
        TTL: int or str
            TTL you want to look at. E.g. "TTL 1", "trial_start" or just "1".
        window: numeric
            Time window both before and after the timestamps
        from_ref: bool
            If you set this to True you'll look at the reference.
        stamps_min: numeric
            Minimum stamp duration to be included
        stamps_max
            Maximum stamp duration to be included
        offset: bool
            Will use the offset of a TTL pulse instead of the onset to
            derive timestamps.

        Returns:
        --------
        A sweepset object with your PE data or a list of these
        objects if TTL refers to multiple channels.
        """
        output = []
        if (TTL.__class__==int) | (TTL.__class__==str):
            TTL = [TTL]

        for ttl in TTL:
            if ttl in self.timestamps.keys():
                temp = self.timestamps[ttl]
                stamps={ttl:{'start':temp}}
            else:
                # This is what we do if we have to derive new stamps
                stamps = self.derive_timestamps(TTL, min_length=stamps_min, max_length=stamps_max,
                        store_stamps = False)
            for key, value in stamps.items():
                if offset:
                    set = self.peri_event(value['end'], window = window, from_ref = from_ref)
                else:
                    set = self.peri_event(value['start'], window = window, from_ref = from_ref)

                # Settings
                set.settings['Z-score'] = True
                set.settings['baseline subtract'] = True
                set.settings['baseline'] = [-10, -2]
                set.settings['display range'] = [-2, window]
                output.append(set)
                set.make_figure()
                plt.suptitle(key, color=self.axcolor)

        # Return
        if len(output) == 1:
            return output[0]
        return output


    def peri_event(self, stamps, window = 10, from_ref = False):
        """
        Will make a peri_event histogram of all data surrounding the stamps.
        
        If you set from_ref = True, it will plot the peri-event traces for
        the reference channel instead, which is a good way to check if
        your effect is due to light or movement artifacts.

        Parameters:
        -----------
        stamps: str or list
            - The timestamps that will be at T=0 in the peri-event histogram. If a string
            will use the timestamps stored in the FIP_signal object by that name.
        window: numerical
            - The window before and after T=0 in seconds.
        from_ref: bool
            - To plot the reference signal instead, set to True

        Returns:
        --------
        A sweepset object containing the PE data

        Examples:
        ---------

        # Make this object, look at the onset of timestamps on 'TTL 1':
        >>> PE = signal.peri_event(signal.timestamps['TTL 1_onset'], from_ref=False)

        # Normalize the PE dataset
        PE.settings['Z-score'] = True

        # Plot the PE data
        PE.make_figure()

        """
        
        # First process the stamps
        if stamps.__class__ == str:
            stamps = self.timestamps[stamps]

        # Make the timeline
        timeline = np.arange(-window, window, 1/self.framerate)

        # Adapt the window
        window = int(window*self.framerate)
          
        # Output variable
        output = np.zeros([len(stamps), len(timeline), self.nr_signals])

        # Original time variable
        original_time = np.zeros([len(stamps), len(timeline), self.nr_signals])
          
        # Grab the data
        if from_ref:
            data = self.get_ref(detrend=True)
        else:
            data = self.get_data()

        # Make sure there is enough data before and after the last stamp
        limit = (window+1)/self.framerate
        indexer = (stamps>data[1, 0, 0]+limit) & (stamps<data[1, -1, 0]-limit)
        if not sum(indexer) == len(stamps):
            print(f'Removing the following timestamps because not enough data is available at the beginning or the end of the recording:')
            print(stamps[~indexer])
            stamps = stamps[indexer]
            output = output[indexer, :, :]
            original_time = original_time[indexer, :, :]

        # Go through all the stamps and grab the correct data
        for i, stamp in enumerate(stamps):
            for signal in range(0, self.nr_signals):
                index = np.argmin(np.abs(data[1,:,signal]-stamp))
                output[i, :, signal] = data[0, index-window:index+window, signal]
                
                # Store the timeline
                original_time[i, :, signal] = data[1, index-window:index+window, signal]
            
        # Make the sweepset object
        new_sweepset = Sweepset(output, timeline, original_time)
        new_sweepset.facecolor = self.facecolor
        new_sweepset.axcolor = self.axcolor

        # Pass on the channel names
        new_sweepset.labels = self.labels        
                      
        return new_sweepset


    def sync_time(self, fip_time:np.array, external_time:np.array):
        """
        Recast the timeline of the FIP signal to align with some external
        timeline.
        
        Note: that sync_time does a simple first-order linear fit so if the 
        discrepancies in timing between the dataset are more complicated, this
        method is not useful. Sync time will keep the frame rate constant 
        throughout the session. 

        Parameters
        ----------
        fip_time : NumPy array
            Set of timestamps taken from the current FIP timeline
        external_time : NumPy array
            Corresponding timestamps to fip_time in the external dataset.

        Returns
        -------
        None.

        """
        
        # Check the input arguments
        if not len(fip_time) == len(external_time):
            print('sync_time only works for equal-length time arrays, please see the docstring')
            return None
        
        # Perform the fit
        p = np.polyfit(fip_time, external_time, 1)
        print(f'Polynomal parameters: {p}')
        
        # Now change the timeline of the signal and the ref
        for i in range(self.signal.shape[2]):
            self.signal[1, :, i] = self.signal[1, :, i]*p[0] + p[1]
            self.raw_ref[1, :, i] = self.raw_ref[1, :, i]*p[0] + p[1]
        
        # The LogAI trace
        self.logAI.loc[:, 'Time'] = self.logAI.Time*p[0] + p[1]

    def _default_formatting(self, ax):
        """
        Just some default formatting to make it look good.
        """

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_facecolor(self.facecolor)
        ax.spines['left'].set_color(self.axcolor)
        ax.spines['bottom'].set_color(self.axcolor)
        ax.xaxis.label.set_color(self.axcolor)
        ax.yaxis.label.set_color(self.axcolor)
        ax.tick_params(colors=self.axcolor)
        ax.title.set_color(self.axcolor)    
                      
        
class Sweepset:
    """
    The sweepset object contains peri-event data. Use the available methods to process and plot.

    Important attributes:
    - RAW_DATA: contains the RAW peri-event data as follows [trials, timeline, signals] however
      it is recommended to use get_data to get the normalized and/or sliced data.
    - SELECTOR: is a list of Bools that can be used to ignore trials. For instance if they contain
      movement artifacts.
    - SETTINGS: contains a dictionary with the settings used to process this peri-event data set.
      for instance, this is where you set the baseline window as well as the Z-score normalization.

    Important methods (see individual doc-strings for more details):
    - GET_DATA: will give you normalized and/or processed data
    - MAKE_FIGURE: will plot the available data
    - GET_AVERAGE: will give you a 1-D average of a channel and a timeline.
    - AUROC: will give you auROC-normalized traces for all channels. This is the preferred method
      of normalizing FIP data if you have enough trials.
    """
    
    def __init__(self, raw_data, timeline, original_timeline):
        # Constructor
        self.raw_data = raw_data
        self.X = timeline
        self.selector = np.ones(raw_data.shape[0], dtype=bool)
        self.nr_signals = raw_data.shape[2]
        self.labels = [f'channel {i+1}' for i in range(raw_data.shape[2])]
        self.original_timeline = original_timeline
        
        # Work on the default settings
        self.settings = {"baseline" : [-5, 0],\
                         "baseline subtract": False,\
                         "Z-score": False,\
                         "display range": [timeline[0], timeline[-1]],\
                         "min-max": False}

        # Some formatting
        self.facecolor = 'w'
        self.axcolor = 'k'
        
    def get_data(self, sliced = False):
        """
        Returns the normalized and baseline-corrected data.

        NOTE: to pick the baseline window and normalization method use the 'settings' property.

        Parameters:
        -----------
        sliced: bool
            - Instead of the returning the complete dataset, will return the part indicated in
            settings['display range'']

        Returns:
        --------
        The output data in a NumPy array like so: [n_trials, n_timepoints, n_channels]

        Examples:
        ---------

            >>> data = get_data() # Get the dataset
            >>> df_data = pd.DataFrame(data[:, :, 0].transpose(),
                                        index = self.X) # Convert channel 1 to a DataFrame
        """
        
        # Grab the raw data
        data = self.raw_data.copy()
        
        # Apply the selector
        data = data[self.selector, :, :]
        
        # Do we do any baseline subtraction
        if (self.settings['baseline subtract']) or (self.settings['Z-score']):
            
            # Figure out the index of the baseline
            b_start = np.argmin(np.abs(self.X - self.settings['baseline'][0]))
            b_end = np.argmin(np.abs(self.X - self.settings['baseline'][1]))

            # Subtract the baseline
            for i in range(0, self.nr_signals):
                    for j in range(0, data.shape[0]):
                        data[j, :, i] = data[j, :, i] - np.mean(data[j, b_start:b_end, i])
                        
                        # Divide by the std if the user requested z-scores
                        if self.settings['Z-score']:
                            data[j,:,i] = data[j, :, i] / np.std(data[j, b_start:b_end, i])
        
        # Should we do any min-max normalization
        if self.settings['min-max']:
            for i in range(0, self.nr_signals):
                for j in range(0, data.shape[0]):
                    # Find the maximum amplitude
                    max_amp = np.max(np.abs(data[j,:,i]))
                    
                    # Divide by the max amplitude
                    data[j,:,i] = data[j,:,i]/max_amp           
        
        # Should we slice the data?
        if sliced:
             # Find out the display range
            range_start = np.argmin(np.abs(self.X - self.settings['display range'][0]))
            range_end = np.argmin(np.abs(self.X - self.settings['display range'][1]))+1
            
            # Slice it
            data = data[:, range_start:range_end, :]
                
        return data
    
    
    def make_figure(self, cmap=None, xlim = None):
        """
        Creates a figure with subplots displaying the data and average traces.

        Parameters:
        -----------
        cmap (matplotlib colormap, optional): The colormap to be used for the heatmap. 
            If not provided, a default colormap is used.
        xlim (tuple, optional): The x-axis limits for the data display range. 
            If not provided, the entire range of the data is displayed.

        Returns:
        --------
        Figure handle

        Note:
        This function relies on the presence of the following attributes in the object:
        - self.X: X-axis values for the data
        - self.nr_signals: Number of signals in the data
        - self.settings: Dictionary containing various settings, including 'display range' and 'Z-score'
        - self.get_data(): Method that retrieves the data to be plotted
        """
        
        # Deal with the default cmap
        if cmap is None:
            cmap = sns.diverging_palette(220, 20, sep=20, as_cmap=True)
        
        # Makes the figure based on the settings
        
        # Grab the data
        data = self.get_data()
        
        # Grab the X axes
        X = np.round(self.X, decimals = 2)
        
        # Find out the display range
        range_start = np.argmin(np.abs(self.X - self.settings['display range'][0]))
        range_end = np.argmin(np.abs(self.X - self.settings['display range'][1]))+1
        
        # Apply the display range to both
        data = data[:,range_start:range_end,:]
        X = X[range_start:range_end]
        
        # Make a figure
        figure, axs = plt.subplots(2, self.nr_signals, tight_layout=False, facecolor=self.facecolor)
        _ = [self._default_formatting(ax) for ax in axs.flatten()]
        figure.set_size_inches([12, 5])

        # This is important to make this work when there is only one signal
        axs = axs.reshape((2, self.nr_signals))
        
        # Plot each signal
        for i in range(0, self.nr_signals):
            
            # Grab the data
            data_df = data[:,:,i]
            data_hm = pd.DataFrame(columns = X, data=data_df)
            
            # Heatmap showing each trial
            temp = sns.heatmap(data_hm, cbar_kws = {'orientation': 'horizontal'},
                        ax=axs[0, i], cmap=cmap, xticklabels=False, center=0)
            axs[0, i].set_ylabel('Trial #')   
            axs[0, i].set_title(self.labels[i])
            temp.collections[0].colorbar.ax.tick_params(axis='x', colors=self.axcolor)
            
            # Grab the data to plot the average
            # (Calculating this yourself is faster then seaborn lineplot)
            mean_data = np.mean(data[:,:,i], axis = 0)
            SEM = np.std(data[:,:,i], axis = 0)/np.sqrt(data.shape[0])
            axs[1, i].plot(X, mean_data)
            
            # Plot the average traces
            axs[1, i].fill_between(X, mean_data+SEM, mean_data-SEM, alpha = 0.2)
            axs[1, i].set_xlabel('Time (s)')
            if self.settings['Z-score']:
                axs[1, i].set_ylabel('Signal (Z-score)')
            else:
                axs[1, i].set_ylabel('Signal')
                
            # Bit of formatting
            ylim = axs[1, i].get_ylim()
            axs[1, i].vlines([0], ymin = ylim[0], ymax=ylim[1],
                             linestyle='--', color='grey')
        return figure          
                
    def get_average(self):
        """
        GET_AVERAGE takes the average over all trials of one channel.

        Parameters:
        -----------
        channel: int
            The channel index (note 0-indexing!) of the channel you are interested in.
        sliced: Bool
            If you set this to True, it will only give you indicated in the "display range" setting

        Returns:
        --------
        A 2D matrix with a timeline in the first column and the mean signal in the 2nd.

        """
        
        # Grab the data
        data = self.get_data()

        # Output
        average = pd.DataFrame(index = self.X, columns = self.labels)

        for channel_i in range(self.nr_signals):
            average.loc[:, self.labels[channel_i]] = np.mean(data[:, :, channel_i], axis = 0)
        
        # Figure out if we should slice
        # Should we slice the data?
        # DO THIS
        # if sliced:
        #      # Find out the display range
        #     range_start = np.argmin(np.abs(self.X - self.settings['display range'][0]))
        #     range_end = np.argmin(np.abs(self.X - self.settings['display range'][1]))+1
            
        #     # Slice it
        #     average = average[range_start:range_end, :]

        return average
    
    
    def get_1d_data(self, sliced = True):
        # Returns a matrix with a the events in a row (instead of stacked on
        # top of each other.)
        
        # Grab the data
        X = self.original_timeline
        Y = self.get_data(sliced = sliced)
        
        # We might have to slice X
        if sliced:
            # Find out the display range
            range_start = np.argmin(np.abs(self.X - self.settings['display range'][0]))
            range_end = np.argmin(np.abs(self.X - self.settings['display range'][1]))+1
            
            # Slice it
            X = X[:, range_start:range_end, :]
        
        # Reshape
        X = X.reshape([-1, self.nr_signals])        
        Y = Y.reshape([-1, self.nr_signals])
        
        return Y, X

    
    def auroc(self):
        """
    	Will normalize the PE_data by calculating the area under the curve of the
    	receiver opperant characteristic at every time point. The final result is
    	a 1-D array with bounds of 0 and 1 and the baseline at 0.5. A vallue of 1
    	at T = X means that there is at least one threshold that perfecty
    	distinguises the values at timepoint X from the values of the baseline. It
    	is a very powerfull method, but only if there are sufficient trials.
    	The input should be a Dataframe with the trials as columns and a timeline 
    	as index. All timepoints with T<0 will be defined as baseline unless a
    	specific baseline is given.
        
    	See Cohen et al. Nature, 2012 for a detailed explanation of this method.
        
            
    	"""
    
        # For every channel        
        # Output data
        column_names = self.labels
        new_data  = pd.DataFrame(index = self.X, columns=column_names)
        
        # Grab all data
        all_data = self.get_data()
        
        # Figure out the baseline
        baseline = self.settings['baseline']
        baseline = (new_data.index>=baseline[0]) & (new_data.index<baseline[1])
        
        for channel_i in range(self.nr_signals):
    
            # Grab the data
            PE_data = pd.DataFrame(index = self.X, 
                   data = all_data[:, :, channel_i].transpose())
            
            # Build a linspace of threshold we will try
            thresholds = np.linspace(PE_data.min().min(), 
                                     PE_data.max().max(), 100)
         
            # Figure out what fraction of baseline above threshold
            base_points = PE_data[baseline].values.reshape(-1)
            l = len(base_points)
            pBaseline = [(sum(base_points>t)/l) for t in thresholds]
      
            # For every timepoint build ROC
            l = PE_data.shape[1]
            for time in new_data.index:
    
                # Pandas indexing is actually really slow, so we want to do it only
                # once every loop.
                temp = PE_data.loc[time, :].values
            
                # Calculate the fraction of trials that the signal is above any threshold
                pSignal = [(temp>t).sum()/l for t in thresholds]
            		
                # Calculate the auROC and store it
                new_data.loc[time, self.labels[channel_i]] =\
                    metrics.auc(pBaseline, pSignal)
            
            # Just making sure it's floats not objects or anything
            new_data = new_data.astype(float)
    
        return new_data
    
    def get_peaks(self, interval):
        
        # Grab the interval
        indexer = (self.X>=interval[0]) & (self.X<interval[1])
        
        # Grab the peaks
        data = self.get_data()
        data = data[:, indexer, :]
        peaks = np.max(data, axis=1)
        
        # columns
        columns = [f'channel_{i}' for i in range(peaks.shape[1])]
        
        return pd.DataFrame(data=peaks, columns = columns)

    def _default_formatting(self, ax):
        """
        Just some default formatting to make it look good.
        """

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_facecolor(self.facecolor)
        ax.spines['left'].set_color(self.axcolor)
        ax.spines['bottom'].set_color(self.axcolor)
        ax.xaxis.label.set_color(self.axcolor)
        ax.yaxis.label.set_color(self.axcolor)
        ax.tick_params(colors=self.axcolor)
        ax.title.set_color(self.axcolor)      
        
                
# This is to run FIP_signal from the command line.       
if __name__ == '__main__':
    import sys

    # EXAMPLE here we use the 2nd argument provided to indicate the dataset.
    if sys.argv[-1][-3:] == 'mat':
        signal = FIP_signal(filename = sys.argv[-1])

        # This is pure vanity
        signal.facecolor = 'k'
        signal.axcolor = 'w'

        plt.ion()
        signal.plot()       