# Fipster
Fipster is a set of Matlab scripts for quick analysis of fiber photometry (FIP) data. To get started simply type:

    >> signal = FIP_signal('User input');

This will open a GUI where one can select a .mat file containing signal (Calcium depended, CD) and reference (405nm emission) traces. There are currently two supported formats. One is that exported by Fipgui (Deisseroth lab) the other is the .mat files exported by the FIP_aquisition class which is also part of this toolbox. If a LogAI file with the same name is found in the same folder, this file will be imported as well.

# Example workflow
Fipster uses a GUI that should be relatively straightforward. The main component, 'FIP_signal' is a Matlab class and supports dot notation for access to non-private properties and methods. This is an example of a standard work flow to obtain time-locked (peri-event) plots of FIP signal.

1. Import signal

		>> signal = FIP_signal('User input');

2. Normalize signal
	Usually the 405nm signal is fitted to the CD signal by adding the 1st and multiplying by the 2th polynomial of the best fit (least squares). This method is called 'polyfit_2' and it is applied by default. Right-click on the reference trace (orange) to select different methods. Additionally, to prevent over-fitting the reference trace is usually smoothed, this is not done by default. By clicking on the 'normalized signal' button of pressing the up or down arrow one can switch to the normalized signal. This is automatically updated. Right-click on the normalized signal to change the signal unit (∆F/F, F or Z-scores). Mostly likely ∆F/F is the preferred signal unit.

3. Apply time offset
	To align different data streams time-stamps are usually collected by the FIP setup using the analog inputs (AI) ports. Both Fipgui and FIP_aquisition store logAI signals in a separate file that ends with 'logAI.csv'. This file is automatically imported by FIP_signal if it is on the same path as the .mat file that contains the FIP data. All collected AI is displayed in the bottom plot (the LogAI plot). Right click on any signal to derive timestamps or (if you are really sure that there is only 1 time stamp in the signal) directly click the option 'use for time alignment'. Otherwise click on the button 'Time stamps' to display time stamps and now click on the time stamp you want to use for time alignment.

4. Import or derive time stamps for peri-event plot.
	Timestamps that are collected elsewhere can be imported as follows:

			>> signal.import_timestamps(stamps,name) 
			stamps: should be 1 x n array of stamps in seconds
			name: should be a char (between ' and ')

	Additionally it is possible to derive timestamps for the LogAI plot (right click) in this case it is possible to skip step 3.

	Note: timestamps are not automatically offset when a new time alignment is applied. So pay attention to the order in which step 3 and 4 are performed. If timestamps from multiple different sources are applied, it might be more efficient to offset them first before integration into FIP signal.

5. Make a peri-event plot
	Right click on the normalized signal and click peri-event plot. Select the time stamps you want to use. This will open a new figure with every trial (sweep) displayed. Note that the event plot is also put in the workspace and that you have access to it's settings by clicking the object 'peri_event'. Additionally most settings are available via a context menu which appears by right-clicking the plot or individual sweeps.

6. Finishing touch
	Remove failed sweeps (e.g. movement artifacts) by clicking on them or pressing 'r' when they are active (switch active sweep using left and right arrow). You can switch between different channels using the space bar. Use the context menu to select a baseline correction method. Note that you should probably not correct baseline if you imported Z-scores (see step 2), but if you have the sweepset object calculate Z-scores, then baseline correction is OK, as it will automatically be applied in the correct order, even it you did not request them in the right order. Note that there are different ways of calculating the z-score, by default the sweepset object uses the baseline window, but you can change this by changing the peri_event.settings.Z_scores_over_baseline boolean to 'false'. 
		
		>> peri_event.settings.Z_scores_over_over_baseline = false;
		
	Finally use 'heatplot' or 'export data' (both in the context menu) to display all selected sweeps. Or export the average trace by clicking 'export average' on the average trace (green).

