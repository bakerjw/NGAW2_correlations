% Call various functions in this directory, to produce the figures in the
% following manuscript:
%
% Baker, J. W., and Bradley, B. A. (2017). "Intensity measure correlations 
% observed in the NGA-West2 database, and dependence of correlations on 
% rupture and site parameters." Earthquake Spectra, 33(1), 145–156.
%
% Note that many of the scripts below start by clearing the workspace,
% so the scripts named below may need to be run individually if you would
% like to analyze intermediate data
%
% Created by Jack Baker
% Last modified June 2, 2016
% Modified December 6, 2021 to update citation information


clear; close all; clc;

%% evaluate GMPE to compute spectral residuals
getResiduals % mixed effects residuals calculation (this is slow)
mergeIMs % merge in non-Sa residuals

%% evaluate correlation predictive models
corrPredictions

%% compare models
cd('compare rho models')
compareRhoModels
cd ..

%% spectral acceleration correlations
corrSa 

%% correlations related to other IMs 
corrNonSa   

%% modify data to enforce positive definiteness
makePosDef

%% M/R/Vs30 dependence
corrVsRup
scatterVsMagDist % scatter plots to see what's driving correlations




