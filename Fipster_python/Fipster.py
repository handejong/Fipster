#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 14:41:27 2019

@author: handejong
"""

# These are neccesary for file browsing
from scipy.io import loadmat
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from scipy import signal
from sklearn import metrics


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
        self.timestamps = []
        self.detrend = False
        self.external_signal = []
        self.smooth = []
        
        # Set the settings.
        self.settings = {
                'fit ref': 'polyfit',
                'signal unit': '∆F/F'}

        
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

                    
    def load_data(self, filename):
        # Load .mat file into the object
        
        print('Loading: {0}'.format(filename))
        
        # Try two different ways of importing old and new Matlab files
        try: 
            data = loadmat(filename)
        except:
            data = h5py.File(filename, mode='r')
              
        # Unpack individual variables  
        for key, value in data.items():
            
            if key == 'framerate':
                self.framerate = np.int(value[()])
            
            if key == 'notes':
                self.notes = value
            
            if key == 'signal':
                self.signal = np.array(value[()])
            
            if key == 'ref':
                self.raw_ref = np.array(value[()])
                
            if key == 'labels':
                self.labels = value
        
        # If there is a ref, divide the framerate by 2
        if not self.raw_ref==[]:
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
            print('Loading: {0}'.format(new_filename))
            self.logAI = pd.read_csv(new_filename)
            
            # Name the columns
            columns = ['Time']
            for i in range(1, self.logAI.shape[1]):
                columns.append('TTL {0}'.format(i))
            self.logAI.columns = columns   
            
            # Remove columns if there is no signal above 1mV
            for i in range(1, len(columns)):
                
                if np.max(self.logAI[columns[i]]) < 1:
                    self.logAI = self.logAI.drop(columns[i], axis = 1)
            
        # Finish up
        self.nr_signals = self.signal.shape[2]
        self.hasdata = True
        self.filename = filename
        self.external_signal = [False] * self.signal.shape[2]
        self.smooth = [False] * self.signal.shape[2]

        
    def import_signal(self, new_signal):
        # Import external signal
        
        # Error handeling
        # TODO
        #   .... does it have it's own timeline?
        #   .... is the sample rate equal to the internal signal
        #   .... does it come with a refference?
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
        # Checks if all settings are vallid (to prevent errors).
        
        # Cycle through the keys and check them all
        for key, value in self.settings.items():
            
            found_key = False
            
            if key == 'fit ref':
                found_key = True;
                if not (value == 'polyfit' or value == 'means' or value =='no fit'):
                    raise Exception("{0} is not a valid fit method. pick from: " \
                                    "\'polyfit\', \'means\', or \'no fit\'.".format(value))
                    
            if key == 'signal unit':
                found_key = True;
                if not (value == '∆F/F' or value == 'Z-score' or value == 'deltaF'):
                    raise Exception("{0} is not a valid fit method. pick from: " \
                                    "\'∆F/F\', \'deltaF\', or \'Z-score\'.".format(value))
                    
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
        
        # Subtract the fitted reference, add the mean vallue
        for i in range(0, self.nr_signals):
            data[0, :, i] = data[0, :, i] - ref[0, :, i] + np.mean(ref[0, :, i])
            
            # should we detrend? (Apply Savitzky-Golay filter)
            if self.detrend:
                data[0, :, i] = data[0, :, i] - \
                signal.savgol_filter(data[0, :, i], 5001, 1) + \
                np.min(data[0, :, i])
                
            # Do we smooth the signal
            if self.smooth[i]:
                data[0, :, i] = signal.savgol_filter(data[0, :, i], 
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
    
    
    def get_ref(self):
        # Calculate the fitted reference signal
        
        # Check settings
        self.check_settings()
        
        ref = self.raw_ref.copy()
        signal = self.signal
        
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
           
            
        return ref
    
    
    def plot(self, signals = 'all', raw_data = True, timestamps = True):
        # Will plot the requested data
        
        # Error handeling on the input arguments
        if len(self.timestamps)==0:
            timestamps = False
        
        # Figure out how many signals to show
        if signals=='all':
            signals = range(0, self.nr_signals)

        # Make a figure
        figure = plt.figure()
        main_plot = plt.subplot(2,1,1)
        
        # Format the figure
        plt.subplots_adjust(top = 0.95, bottom = 0.1, left = 0.1, right = 0.98, hspace = 0.3)
        
        # Plot the signals
        if raw_data:
            ref = self.get_ref()
            for i in signals:
                plt.plot(self.signal[1, :, i], self.signal[0, :, i], lw = 0.1)
                plt.plot(ref[1, :, i], ref[0, :, i], lw = 0.1)
            
            # Format the plot
            plt.xlabel('Time (s)')
            plt.ylabel('Signal (F)')
        
        # Plot normalized data instead
        else:
            data = self.get_data()
            for i in signals:
                plt.plot(data[1, :, i], data[0, :, i], lw = 0.1)
                
            # Format
            plt.xlabel('Time (s)')
            plt.ylabel(self.settings['signal unit'])
        
        plt.title(self.filename)
            
        # Plot the logAI data
        logAI_plot = plt.subplot(2,1,2, sharex=main_plot)
        if not timestamps:
            for i in range(1, self.logAI.shape[1]):
                plt.plot(self.logAI.iloc[:, 0], self.logAI.iloc[:, i],\
                         label = self.logAI.columns[i])
            plt.ylabel('V (mV)')    
            plt.xlabel('Time (s)')
            plt.legend(loc = 'upper right')
        else:
            for i in range(0, len(self.timestamps)):
                plt.scatter(self.timestamps[i]['start'],\
                            np.ones(self.timestamps[i]['start'].shape)*i,\
                            label = self.timestamps[i]['name'])
            plt.ylabel('Stamp #')    
            plt.xlabel('Time (s)')
            plt.legend()
        
        # Return the figure handle
        return figure
    
    
    def derive_timestamps(self, TTL = 'all', names = 'default', min_length = 0,\
                          max_length = 9999):
        # Derive timestamps from the logAI signals
        
        # Output variable
        output_stamps = {}
        

        # Individual columns or all columns
        if TTL == 'all':
            columns = self.logAI.columns[1:]
        else:
            if TTL.__class__ == str:
                columns = [TTL]
            else:
                columns = TTL    
            
        # Figure out TTL stamps names
        if names=='default':
            names = columns
            
        # Just like with columns, make sure it's a list
        if names.__class__ == str:
            names = [names]
            
        # Error handeling
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
            threshold = 1
            TTL_high = signal>threshold
            dif = np.append([0], np.diff(TTL_high.astype(int)))
            
            # Grab the starts and ends
            starts = timeline[dif==1]
            ends = timeline[dif==-1]
            
            # Error handeling
            if ends[0]<starts[0]:
                ends = ends[1:]
            if starts[-1] > ends[-1]:
                starts = starts[:-1]
            if not len(starts) == len(ends):
                print('Unequal number of TTL starts and stops, this is an error!')
            
            # Make the set of stamps
            stamps['start'] = starts
            stamps['end'] = ends
            stamps['TTL'] = col
            stamps['name'] = name
            self.timestamps.append(stamps)
            print('Importing {0} stamps, named: {1}'.format(len(stamps['start']), stamps['name']))
            
            # Can we calculate the stamp length?
            stamps['length'] = stamps['end'] - stamps['start']
            selector = stamps['length']>min_length
            selector[stamps['length']>max_length] = False
            
            # Apply selector
            stamps['start'] = stamps['start'][selector]
            stamps['end'] = stamps['end'][selector]
            stamps['length'] = stamps['length'][selector]

            # Add these stamps to the output
            output_stamps[name] = pd.DataFrame(stamps)
            
        return output_stamps
    
    
    def peri_event(self, stamps, window = 10):
        # Will make a peri_event histrogram of all data surounding the stamps
        # The indexer is a bool that indexes all the actual data used
          
        # Make the timeline
        timeline = np.arange(-window, window, 1/self.framerate)
        
        # Adapt the window
        window = np.int(window*self.framerate)
          
        # Output variable
        output = np.zeros([len(stamps), len(timeline), self.nr_signals])
        
        # Original time variable
        original_time = np.zeros([len(stamps), len(timeline), self.nr_signals])
          
        # Grab the data
        data = self.get_data()
        
        # The indexer
        #indexer = np.zeros([data.shape[1], data.shape[2]]).astype(bool)
        
        # Is there enough data at the end of the stamp?
        data_cutoff = np.max(stamps) + window
        if data[1, -1, 0]<data_cutoff:
            print('Deleting last stamp, because not enough data available.')
            stamps = stamps[:-1]
            output = output[:-1, :, :]
        
        # Go through all the stamps and grab the correct data
        for stamp in range(0, len(stamps)):
            for signal in range(0, self.nr_signals):
                index = np.argmin(np.abs(data[1,:,signal]-stamps[stamp]))
                output[stamp, :, signal] = \
                    data[0, index-window:index+window, signal]
                
                # Store the original timepoints
                original_time[stamp, :, signal] = \
                    data[1, index-window:index+window, signal]
                
                # Select that same data in the indexer
                #indexer[index-window:index+window, signal] = True
            
        # Make the sweepset object
        new_sweepset = Sweepset(output, timeline, original_time)        
                      
        return new_sweepset
                      
        
class Sweepset:
    # ....
    
    def __init__(self, raw_data, timeline, original_timeline):
        # Constructor
        self.raw_data = raw_data
        self.X = timeline
        self.selector = np.ones(raw_data.shape[0], dtype=bool)
        self.nr_signals = raw_data.shape[2]
        self.original_timeline = original_timeline
        
        # Work on the default settings
        self.settings = {"baseline" : [-5, 0],\
                         "baseline subtract": False,\
                         "Z-score": False,\
                         "display range": [timeline[0], timeline[-1]],\
                         "min-max": False}
        
    def get_data(self, sliced = False):
        # Grabs the data based on the raw_data and the settings
        
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
    
    
    def make_figure(self, cmap=None):
        
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
        figure, axs = plt.subplots(2, self.nr_signals, tight_layout=False)
        figure.set_size_inches([12, 5])
        
        # Plot each signal
        for i in range(0, self.nr_signals):
            
            # Grab the data
            data_df = data[:,:,i]
            data_hm = pd.DataFrame(columns = X, data=data_df)
            
            # NOTE: dataFrame for the heatmap so we can do the correct
            # x-axis scaling in a future version.
            
            # If using the default
            # Heatmap showing each trial
            sns.heatmap(data_hm, cbar_kws = {'orientation': 'horizontal'},
                        ax=axs[0, i], cmap=cmap, xticklabels=False, center=0)
            axs[0, i].set_ylabel('Trial #')   
            #axs[0, i].set_xticks(X)
            axs[0, i].set_title('Channel {0}'.format(i+1))
            
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
                
        # Force plot
        plt.show()
                
                
    def get_average(self, channel, sliced = False):
        
        # Grab the data
        data = self.get_data()
        
        # Grab the Timeline
        X = np.expand_dims(np.round(self.X, decimals = 2), axis = 1)
        
        # Make the average (including it's timeline)
        average = np.expand_dims(np.mean(data[:, :, channel], axis = 0), axis = 1)
        average = np.concatenate((X, average), axis = 1)
        
        # Figure out if we should slice
        # Should we slice the data?
        if sliced:
             # Find out the display range
            range_start = np.argmin(np.abs(self.X - self.settings['display range'][0]))
            range_end = np.argmin(np.abs(self.X - self.settings['display range'][1]))+1
            
            # Slice it
            average = average[range_start:range_end, :]
            
        
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
        column_names = [f'channel_{i+1}' for i in range(self.nr_signals)]
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
                new_data.loc[time, f'channel_{channel_i+1}'] =\
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
        
        
                
        
        