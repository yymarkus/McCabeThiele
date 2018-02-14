%Created by: Yonatan Markus
%January 2018
%Purpose: 
%(1) Creating McCabe-Thiele diagrams for distillation columns.
%(2) Predicting NTU (Number transfer units) for column, as well as the
%transition point between stripping and enriching sections (xtrans_new)

%Necessary files:
%PackedDistTransition - Main script for calculations and plotting
%ControlFile - Iterating through PackedDistTransition multiple times for
%different experimental column conditions
%yEnriching - Equation for enriching operating line
%yStripping - Equation for stripping operating line
%fstar - Equation for ethanol-water equilibrium curve
%xb - Predicts bottoms product mol% based on feed and distillate mole flowrates
%and mol% compositions
%NTUenr - Calculates NTU for enriching section
%NTUstr - Calculates NTU for stripping section

%If xb (bottoms mol%) is unknown, the file xb.m should be used to for this
%prediction. To use this calculated xb value, three things must be changed:
%(1) In PackedDistTransition.m, the line '%xb_1 = xb(xf, xd, F, D);' should
%be uncommented.
%(2) In yStripping.m, the y should be calculated with the equation that does
%NOT include xb
%(3) In NTUstr.m, the denominator should be calculated with the equation
%that does not include xb

%% Options
%Saving the output figures to your drive
Saving = false;
SaveData_to_CSV = false;
NTU = false;
Import_data = false;
CSV_filename = '/Users/Yonatan/Documents/University/Year3/2018W/CHE305AppliedChemistryLabIV/Lab4_PackedDistillation/Distillationfile.csv';

%% Import the data
if Import_data == true
    filename = '/Users/Yonatan/Documents/University/Year3/2018W/CHE305AppliedChemistryLabIV/Lab4_PackedDistillation/DistillationRawData.xlsx';
    first_row = 41;
    first_column = 13;
    [~, ~, raw] = xlsread(filename,'Sheet1');
    raw = raw(first_row:end,first_column:first_column+6);

%% Create output variable
    data = reshape([raw{:}],size(raw));

%% Allocate imported array to column variable names
    R0 = data(:,1);
    xd0 = data(:,2);
    xf0 = data(:,3);
    xb0 = data(:,4);
    D0 = data(:,5);
    F0 = data(:,6);
    B0 = data(:,7);

%% Clear temporary variables
    clearvars data raw;
end
%% Manually Input Variables 
 R0 = [1; 1; 1; 1; 2; 3];
 xd0 = [0.7550182623; 0.7082146251; 0.6828964658; 0.4211082015; 0.2127421787; 0.5987939856];
 xb0 = [0.005962989282; 0.005314175964; 0.006179549448; 0.002505749081; 0.001250434856; 0.01136074851];
 xf0 = [0.08615848586; 0.0713588374; 0.09064226432; 0.07241758752; 0.07545118289; 0.08792153704];
 F0 = [0.004625271467; 0.00948343042; 0.00475190704; 0.009807787257; 0.009747225206; 0.009559981999];
 D0 = [0.001000524914; 0.001039757953; 0.002028966601; 0.002104185912; 0.00186101202; 0.000924866794];
 B0 = [0; 0; 0; 0; 0; 0.1];
q0 = [1.5; 1.5; 1; 1; 1; 1];

xtrans = zeros(length(xf0), 1);
NTU = zeros(length(xf0), 1);

for i = 1:length(xf0)
     [xtrans(i) NTU(i)] = PackedDistTransition(xf0(i), xd0(i), xb0(i), R0(i), F0(i), D0(i), B0(i), i, Saving, NTU, q0(i));
end

%% HTU
z = 1.72; %m
HTU = z./NTU;

%% Saving
if SaveData_to_CSV == true
    dlmwrite(CSV_filename,['xtrans'], 'delimiter', ' ', '-append');
    dlmwrite(CSV_filename,['NTU'], 'delimiter', ' ', '-append', 'coffset', 1);
    dlmwrite(CSV_filename,[xtrans NTU], '-append');
end
    
    