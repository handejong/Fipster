# Fipster
Fipster is a set op Matlab scripts for quick analysis of fiber photometry (FIP) data. To get started simply type:

    >> signal = FIP_signal('User input');

This will open a GUI where one can select a .mat file containing signal (Calcium dependend) and reference (405nm emmission) traces. The only currently suported format is those exported by Fipgui (Deisseroth lab). If a LogAI file with the same name is found in the same folder, this file will be imported as well.

# Example workflow
Fipster uses a GUI that should be relatively straighforward. 
