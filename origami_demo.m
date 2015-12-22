function origami_demo(type)

% This DEMO illustrates how to use the ORIGAmI package.
% INPUT: 'type' denotes a string argument to select which program to demo:
% 'tex' = demo getTEX86.m
% 'alk' = demo getAlkenones.m
% 'gc' = demo getGCAreas.m
% 'hplc' = demo getHPLCArea.m

if strcmpi(type,'tex')
    display('Calculating TEX86 and BIT from a sample without an internal standard')
    %define retention times
    RT={21.98;21.01;20;[15.62 16.91];14.14;13;11.83;10.88};
    TEX_example1=getTEX86('Folders','TEX_nostd','standard',0,'rt',RT,'showgraphs',1)
    
    display('Calculating TEX86 and BIT from a sample with an internal standard')
    %this time use default RTs in code, which are close to the sample RTs.
    %note that standard peak area appears last.
    TEX_example2=getTEX86('Folders','TEX_std','standard',1,'showgraphs',1)
    
elseif strcmpi(type,'alk')
    display('Calculating UK indices from a GC trace, with an internal standard at 33.6 minutes')
    %define retention times.
    RT=[33.66,36.66,37.25,39.8,40.06,40.48,40.66];
    UK_example1=getAlkenones('files','Alkenone_1.csv','rt',RT,'showgraphs',1)

elseif strcmpi(type,'gc')
    %integrate the even homologs of a FAME trace.
    display('Calculating GC peak areas from a trace of Fatty Acid Methyl Esters')
    %define retention times.
    RT=[13.5,15.15,16.9,18.64,20.35,21.99,23.55,25.01,26.46];
    %note smoothparam is set to a small number for this file.
    GC_example1=getGCAreas('files','FAMEs.csv','rt',RT,'showgraphs',1,'smoothparam',1)
    
elseif strcmpi(type,'hplc')
    %integrate a large set of isoGDGTs and brGDGTs.
    display('Calculating HPLC peak areas from multiple EICs')
    %define files to integrate
    files={'744.txt';'1302.txt';'1300.txt';'1298.txt';'1296.txt';'1292.txt';'1050.txt';'1048.txt';'1036.txt';'1034.txt';'1032.txt';'1022.txt';'1020.txt';'1018.txt'};
    %define RTs, including some multiple peaks for certain files.
    RT={29.6;16.3;17.8;19.7;21.5;[23.65,25.4];[38.3,39.1];[39.4,40.4];[40.9,41.3];[42,42.5];[43.1,43.6];43.4;44.5;45.6};
    %note window setting is small to ensure accuracy and speed.
    HPLC_example1=getHPLCAreas('folders','HPLC_data','files',files,'rt',RT,'window',2,'showgraphs',1)
else
    error('%s is not a recognized demo type',type)
end