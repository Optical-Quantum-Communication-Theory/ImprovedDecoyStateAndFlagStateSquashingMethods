%% FUNCTION NAME: Main
% Main entry point function for key rate calculation.
% Please find manual for detailed explanations on the functions and the software structure.
%
% The user can start the program by choosing a preset to run or modifying a custum preset.
% Each preset contains three functions 
% 'setDescription', 'setParameters', and 'setOptions', 
% where respectively the protocol/channel description files, protocol parameters, 
% and software/solver options are inputted.
%%

format long
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%% Setting User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

preset='pmSixStateWCP_decoy';
[protocolDescription,channelModel,leakageEC,parameters,solverOptions]=feval(preset);

%%%%%%%%%%%%%%%%%%%%% Run Main Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call main iteration function
results=mainIteration(protocolDescription,channelModel,leakageEC,parameters,solverOptions);

%%%%%%%%%%%%%%%%%%%%% Output Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save the results to file
%save('Outputs/output_6state_reduced2.mat','results','parameters');
%save('output.mat','results','parameters');

% %can also load a previous session's result to plot it
% %(can comment out main iteration above to skip computation)
% load('output.mat','results','parameters_scan')

%can uncomment this line to output debugging info
% results.debugInfo

%automatically parse and plot the results (optional)
%the third optional argument is the plotting style
%available options for 1D data:
%1.'linear': plot x and y linearly
%2.'linear-log': plot x versus log10(y)
%3.'km-log': plot -log10(x)*10/0.2 (x is assumed to be transmittance and converted to km) versus log(y)
%4.'none': do not plot
plotResults(results,parameters,'linear-log')