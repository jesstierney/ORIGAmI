# ORIGAmI
A Matlab package for automated peak integration of GC and LC data.

ORIGAmI
This ReadMe covers the usage of the main functions in the ORIGAmI (ORganIc Geochemistry peAk Integration) package.

Installation
Put the ORIGAmI folder in your main Matlab directory (e.g., Documents/MATLAB) and update your path in Matlab to include the folder and all its functions.

getTEX86

getTEX86 calculates the peak areas, TEX86 and BIT indices from HPLC traces. We assume that each sample has its own folder containing all relevant EIC files as .txt files. E.g., the folder ‘Sample1’ would contain a text file for each ion (‘1302.txt’,’1300.txt’, etc…). If you use Chemstation, you can use the macro txtfile.mac included in this package to export Chemstation files into the correct format.
The inputs and outputs of getTEX86.m are:

Output=getTEX86(varargin)

INPUTS: varargin - specifies optional parameter name/value pairs. ‘Folders’ - names of folders to calculate indices for. For one folder, this input can be a single string. To calculate the areas for multiple files at once, use a cell array of strings, e.g. {‘Sample1’;’Sample2’} Default: Integrates all folders in the current directory. ‘Files’ - names of individual files to calculate areas for. The files should be a cell array of strings, Default: {‘1022.txt’,’1036.txt’,’1050.txt’,’1292.txt’,’1296.txt’,’1298.txt’,’1300.txt’,’1302.txt’,’744.txt’}; * Note: 744 is the EIC corresponding to a commonly-used internal standard (see Huguet et al., 2006). To include the standard by default, you must set the value for 'Standard' to 1. ‘RT’ - Estimated retention times of the peaks of interest in each file, entered as a cell array. The order should be the same as files, e.g.: RT={1022,1036,1050,1292 (two RTs),1296,1298,1300,1302,744}; Example: RT={21.8;20.9;19.9;[15.7 17.1];14.2;13.06;11.86;10.91;17.02}; Note: You can change also change RTs in the code to match your mass spec values by changing the variable 'RT_Sample' in the subfunction 'getPeaks’. The code uses these values as the default. ‘Standard’ - 0 or 1. Choose 1 to include the 744 standard in the default files, and 0 to not include 744. Default : 0 ‘ShowGraphs’ - 0 or 1, where 1 produces a set of graphs for each sample and each EIC that show the EIC trace and the Gaussian fit. Default: 0 ‘Window’ - Peaks must fall within a window of length w that is centered around the retention times defined by RT. To improve the speed of the calculation, make the window smaller. To improve the accuracy, you may need to increase the window size. Default: 5

OUTPUT A structure with the following fields: .Folders - List of the folders included, as input .Files - List of the files included, as input .InputTimes - List of the retention times, as input .PeakTimes - Retention times of each peak for each folder .Areas - Areas of each peak for each folder .TEX86 - TEX86 indices corresponding to each folder .BIT - BIT indices corresponding to each folder

NOTES: 1. When calculating the TEX86 index, we assume that we are using the first peak for 1298, the first peak for 1296, the second peak for 1292, and the first peak for 1300. When calculating the BIT index, we assume that we are using the first peak in the files for 1050, 1036, 1022, and 1291. If this is incorrect, you will need to calculate the TEX86 and BIT indices manually with the correct peaks. 2. You will either need to change the default retention times in the code to match typical values from your mass spec (to do this, edit the variable ‘RT_Sample’ with values from one of your samples) or enter the retention times manually using ‘RT’ (recommended).

EXAMPLES
	1.	Calculate areas, TEX86 and BIT for all folders in the current directory, using all default values. Output=getTEX86; 
	2.	Calculate areas, TEX86 and BIT for all folders in the current directory, using default values for files and retention times, but also including the 744 standard and showing the resulting graphs. Output=getTEX86('Standard',1,'ShowGraphs',1); 
	3.	Calculate areas, TEX86 and BIT for two folders, using the same retention times for both folders. Output=getTEX86('Folders',{'ARC_1';'ARC_2'},'RT',{21.54; 20.5; 19.44; [15.29,16.60]; 13.84; 12.72; 11.58; 10.66}); 
	4.	Calculate areas, TEX86 and BIT for two folders, using the different retention times for each folders. Output=getTEX86('Folders',{'ARC_1';'OCE205_1'},'RT’,{{21.54; 20.5; 19.44; [15.29,16.60]; 13.84; 12.72; 11.58; 10.66} ; {21.98;21.006;20.002;[15.623,16.934];14.144;13.001;11.829;10.881}}); 
	5.	Calculate areas, TEX86 and BIT for all folders, using the same retention times for all folders, and an internal standard: Output=getTEX86('RT',{21.8;20.9;19.9;[15.7 17.1];14.2;13.06;11.86;10.91;17.02},'Standard',1); 

getAlkenones

getAlkenones calculates the peak areas, UK37, UK38_et, and UK38_me indices from GC traces. We assume that each sample is a .csv file with time in the first column and intensity in the second column.

The inputs and outputs of getAlkenones.m are:

Output=getAlkenones(varargin)

INPUTS: varargin specifies optional parameter name/value pairs. ‘Files’ - names of files to calculate indices for. If there is only one file, the input can be a string. If there are multiple files, the input should be a cell array of strings. Default: All files in the current directory e.g. {'GC27_1.csv','GC27_3.csv','GC27_4.csv',... 'KNR197_1.csv','KNR197_2.csv','KNR197_3.csv',... 'P178_1.csv','P178_2.csv','P178_3.csv'}; ‘RT’ - Estimated retention times of the peaks of interest in each file, entered as a cell array or double. The order should correspond to the order of elution, e.g. RT={C37:3;C37:2;C38:3et;C38:3me;C38:2et;C38:2me}; or, if you have an internal std: RT={IS;C37:3;C37:2;C38:3et;C38:3me;C38:2et;C38:2me}; Example: RT = [33.5,37.7,38.4,41.0,41.2,41.7,41.9]; You may input 1-2 internal standards, as long as they elute before the alkenones. Note: You can change also change RTs in the code to match your GC values by changing the variable 'RT_Sample' in the subfunction 'getPeaks’. The code uses these values as the default. ‘ShowGraphs’ - 0 or 1, where 1 produces a graph for each sample that shows the trace, smoothed trace, and the Gaussian fits. Default: 0 ‘Window’ - Peaks must fall within a window of length w that is centered around the retention times defined by RT. To improve the speed of the calculation, make the window smaller. To improve the accuracy, you may need to increase the window size. Default: 2 ‘SmoothParam’ - Amount to smooth the data before integration. Use lower values (e.g., 2) for larger peaks to avoid a reduction in peak height. Default: 15

OUTPUT A structure with the following fields: .Files - List of the files included, as input .InputTimes - List of the retention times, as input .PeakTimes - Retention times of each peak for each file .Areas - Areas of each peak for each file .UK37 - UK37 indices corresponding to each file .UK38_et - UK38_et indices corresponding to each file .UK38_me - UK38_me indices corresponding to each file

NOTES:
	1.	When calculating the UK37, UK38_et, and UK38_me indices, we assume that the peaks used to calculate these indices are the final six peaks in each file. If this is incorrect, you will need to calculate the indices manually with the correct peaks.
	2.	You will either need to change the default retention times in the code to match typical values from your GC (to do this, edit the variable ‘RT_Sample’ with values from one of your samples) or enter the retention times manually using ‘RT’ (recommended).

EXAMPLES
	1.	Calculate areas, UK37, UK38_et, and UK38_me for all files in the current directory, using all default values. Output=getAlkenones; 
	2.	Calculate areas, UK37, UK38_et, and UK38_me for one file, using default values for retention times, but showing the resulting graph. Output=getAlkenones('Files','GC27_2.csv','ShowGraphs',1); 
	3.	Calculate areas, UK37, UK38_et, and UK38_me for all files in the current directory, using the same retention times for all folders, and showing the resulting graphs. Output=getAlkenones('RT',[36.62, 37.25, 39.78, 40.05, 40.47, 40.65],'ShowGraphs',1); 
	4.	Calculate areas, UK37, UK38_et, and UK38_me for two files, using the different retention times for each folders, and showing the resulting graphs. Output=getAlkenones('Files',{'GC27_1.csv','KNR197_1.csv'},'RT',{[36.62, 37.25, 39.78, 40.05, 40.47, 40.65];[33.63, 36.66, 37.19, 39.78, 40.05, 40.42, 40.62]},'ShowGraphs',1); 

getHPLCAreas

getHPLCAreas calculates the peak areas of high-performance liquid chromatography data. We assume that each sample has its own folder containing all relevant EIC files. We assume that the files are .txt files. If you use Chemstation, use the txtfile.mac macro in this package to export data to the correct format.
The inputs and outputs of getHPLCAreas.m are:

Output=getHPLCAreas(varargin)

INPUTS: varargin - specifies optional parameter name/value pairs. ‘Folders’ - names of folders to calculate indices for. For one folder, this input can be a single string. To calculate the areas for multiple files at once, use a cell array of strings, e.g. {‘Sample1;’Sample2’} Default: All folders in the current directory ‘Files’ - names of individual files to calculate areas for. The files should be a cell array of strings, e.g. {'654.txt';'744.txt';'1302.txt';'1300.txt';... '1298.txt';'1296.txt';'1292.txt';'1050.txt';... '1048.txt';'1036.txt';'1034.txt';'1022.txt';... '1020.txt';'1018.txt’}; Default: All files in the folder that have a retention time defined in the variable 'RT_Sample' ‘RT’ - Estimated retention times of the peaks of interest in each file. The number and order should match ‘Files’. Default: {16.85;21.54;20.5;19.44;[15.29,16.60];13.84;12.72;11.58;10.66}, corresponding to the files: {'744';'1022';'1036';'1050';'1292';'1296';'1298';'1300';'1302'} Note: You can change this default value to match your mass spec values by changing the variable 'RT_Sample' ‘ShowGraphs’ - 0 or 1, where 1 produces a set of graphs for each sample and each EIC that show the EIC trace and the Gaussian fits. Default: 0 ‘Window’ - Peaks must fall within a window of length w that is centered around the retention times defined by RT. To improve the speed of the calculation, make the window smaller. To improve the accuracy, you may need to increase the window size. Default: 5

OUTPUT A structure with the following fields: .Folders - List of the folders included, as input .Files - List of the files included, as input .InputTimes - List of the retention times, as input .PeakTimes - Retention times of each peak for each folder .Areas - Areas of each peak for each folder

NOTES: 1. You can change RT_Sample in the code to suit your needs, but we recommend entering a cell array of ‘Files’ and a matching cell array of ‘RT’ every time you use the function, for ease. This function can integrate multiple peaks per EIC, if needed.

EXAMPLES
Note that usage is similar to getTEX86.
1.	Calculate areas for all files* in all folders in the current directory, using pre-determined files and RT variables, and a custom window size:

output=getHPLCAreas('files',files,'rt',RT,'window’,2);

getGCAreas

getGCAreas calculates the peak areas of gas chromatography data. We assume that each sample is a .csv file with time in the first column and intensity in the second column.

The inputs and outputs of getGCAreas.m are:

Output=getGCAreas(varargin)

INPUTS: varargin specifies optional parameter name/value pairs. ‘Files’ - names of files to calculate indices for. If there is only one file, the input can be a string. If there are multiple files, the input should be a cell array of strings. Default: All files in the current directory e.g. {'GC27_1.csv','GC27_3.csv','GC27_4.csv',... 'KNR197_1.csv','KNR197_2.csv','KNR197_3.csv',... 'P178_1.csv','P178_2.csv','P178_3.csv'}; ‘RT’ - Estimated retention times of the peaks of interest in each file. Default is something like alkenone peaks: [39.75,40.38,43.08,43.34,43.75,43.96] We recommend plotting a test sample and picking out the RTs of the peaks you want to integrate, then saving those values as a double for input. ‘ShowGraphs’ - 0 or 1, where 1 produces a graph for each sample that shows the trace, smoothed trace, and the Gaussian fits. Default: 0 ‘Window’ - Peaks must fall within a window of length w that is centered around the retention times defined by RT. To improve the speed of the calculation, make the window smaller. To improve the accuracy, you may need to increase the window size. Default: 2 ‘SmoothParam’ - Amount to smooth the data before integration. Use lower values (e.g., 2) for larger peaks to avoid a reduction in peak height. Default: 15

OUTPUT A structure with the following fields: .Files - List of the files included, as input .InputTimes - List of the retention times, as input .PeakTimes - Retention times of each peak for each file .Areas - Areas of each peak for each file

EXAMPLES
Note that usage is similar to getAlkenones.
	1.	Calculate areas for all files in all folders in the current directory, using pre-determined RT variables, and a custom window size:

output=getGCAreas('rt',RT,'window’,2);

showplot

This simple function allows the user to quickly plot GC and HPLC traces (.csv or .txt format).

EXAMPLES
Plot a whole GC trace: showplot(‘Sample1.csv’);
Plot a HPLC EIC: showplot(’1292.txt’);
Plot part of a GC trace from 10 to 30 minutes: showplot(‘Sample1.csv’,[10 30]);
