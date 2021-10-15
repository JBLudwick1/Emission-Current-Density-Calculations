%This is the main script to call the lookup table generator
%--------------------------------------------------------------------------
% Input parameters must be in the same array order for UBs and LBs
%  LBs: lower bounds for [input parameter 1, input parameter 2, ... , input parameter N]
%  UBs: upper bounds for [input parameter 1, input parameter 2, ... , input parameter N]
%  Array parameter order must be the same as numerical function input parameter order
%  If a parameter is stay constant during the sweep, set the upper and
%  lower bounds to be equal
%--------------------------------------------------------------------------
% Cell arrays of input parameters to vary in a special function space
% 0 value indicates to treat quantity in linear space
% These cell arrays must be in the same order as UBs/LBs
% special_funcs_inv should be the inverse funcitons of special_funcs
% These special functions should be able to handle 1D array inputs
% This is useful with parameters that behave much less linearly in certain
% regions. For example, radius of curvature is far less linear for small
% values, and therefore can more efficently be used in log10 space.
%--------------------------------------------------------------------------
% Conduct linear interpolation in log space? Quantities such as emission
% current density tend to perform better in log space.
% 0: Coduct standard linear interpolation
% 1: Coduct linear interpolation in log space
%--------------------------------------------------------------------------
%Numerical function to evaluate (name of matlab funciton in string format)
%This function must allow 1D array inputs for each input parameter
%This function must output an N-Dimensional matrix where N is the quantity
%of input parameters. All combinations of the input parameters are to be
%calculated.

% Input parameters
LBs = [0.5,300,20,2,5]; %Lower Paramater Bounds. This example contains the parameters, in order from left to right, electric field magnitude [V/nm], temperature [K], radius of curvature [nm], work function [eV], fermi energy [eV]
UBs = [10,3000,5000,6.5,5]; %Upper Paramater Bounds. This example contains the same parameters as LBs

% Cell arrays of input parameters to vary in a special function space
special_funcs = {0,0,@(x)log10(x),0,0};
special_funcs_inv = {0,0,@(x)10.^x,0,0};

%Max relative error threshold across entire lookup table
global error_threshold;
error_threshold = 10; %This example sets a 10% relative error maximum

% Conduct linear interpolation in log space?
logspace_interpolation = true;

output_full_path = ''; %Full path (including file name) for finished lookup table

%Numerical function to evaluate (name of matlab funciton in string format)
num_fun_str = 'NUMERICAL_FUNC.QMWI'; %This example uses the QMWI function

%Call lookup table generator
LUT_gen.lookup_table_generator(LBs,UBs,special_funcs,special_funcs_inv,logspace_interpolation,output_full_path,num_fun_str)