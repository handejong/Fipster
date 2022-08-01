#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 19:21:53 2021

@author: han
"""

# Import libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Custom libraries
import Fipster as fip

# Settings
plt.ion()

# <codecell> run analysis on one mouse
dataset = '../raw_data/example_1.mat'

# Load data
signal = fip.FIP_signal(filename = dataset)

# Export the raw data to CSV for analysis with 
raw_data = signal.get_data()
raw_data = pd.DataFrame(index = raw_data[1, :, 0],
    columns = ['Channel_1', 'Channel_2', 'Channel_3'],
    data = raw_data[0, :, :])
raw_data.index.name='Time (s)'
raw_data.to_csv('example_data.csv')
print(f'Raw data exported to .csv file.')

# Have a look at the signal
signal.plot(raw_data=False)
plt.show()

# Convert the TTL signal (bottom pannel) to timestamps
licks = signal.derive_timestamps('TTL 1',
    names=['Licks'])['Licks']
cue_onset = signal.derive_timestamps('TTL 2', 
    min_length=0.9, 
    max_length=1.1,
    names=['CS onset'])['CS onset']

# Make some peri_event plots
def make_peri_event(stamps):
    
    # Make the sweeset object
    sweepset = signal.peri_event(stamps)
    
    # Include some settings
    sweepset.settings["baseline"] = [-10, 0]
    sweepset.settings["Z-score"] = True
    
    return sweepset

# Make and plot a peri-event dataset
CS_onset = make_peri_event(cue_onset.end)
CS_onset.make_figure(); plt.suptitle('CS onset')

# <codecell> extract data for further analysis
def get_data(sweepset):
    
    # hardcoded fibers
    fibers = ['Channel_1', 'Channel_2', 'Channel_3']
    
    # Output
    output = pd.DataFrame(columns = fibers)
    for i, fiber in enumerate(fibers):
        output.loc[:, fiber] = sweepset.get_average(i)[:, 1]
        
    output.index = sweepset.get_average(i)[:, 0]
    
    return output

# Export the peri-event data
response_to_CS = get_data(CS_onset)

# Quick inspection of the dataframes
response_to_CS.plot(title='Response to CS', ylabel='Z-score',
                    xlabel='Time (s)')

# auROC plot
auROC_data = CS_onset.auroc()
auROC_data.plot() 
plt.title('auROC example')