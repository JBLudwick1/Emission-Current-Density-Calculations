clear; clc; close all;

%--------------USER INPUT--------------
%Enter upper (UB) and lower (LB) bounds for each input parameter
%input parameters must be in the same array order for UBs and LBs
LBs = [1e6,300]; %lower bounds for [input parameter 1, input parameter 2, ... , input parameter N]
UBs = [2e8,3000]; %upper bounds for [input parameter 1, input parameter 2, ... , input parameter N]
%Above example is for a 2D lookupt table with an Electric Field range of 1e6 to 2e8 V/cm and 300 to
%3000 K

%For percent error, a value of 1 here would represent a max error of 1%
error_threshold = 1;

%error_type 0: percent error
%error_type 1: absolute error
error_type = 0;

%Enable (1) [Disable (0)] logspace_interpolation to measure errors by linearly interpolating a function in
%logspace. This can reduce lookup table sizes depending on the function.
%The resulting table must always use linear inerpolation in logspace in
%subsequent simulations if this value is enabled.
%For emission current density, it is beneficial to enable logspace_interpolation
logspace_interpolation = 1;

%Filepath to store generated table
output_table_savepath = '';

%Table file name
fn = 'Table';

%FINAL USER INPUT INFO-----------------------------------------------------
%Lastly, the user must specify the function to numerically solve in the
%'function_calculation' function at the bottom of this code. This current
%example calls a numerical determination of emission current density.

%Error estimates greater than 100% error may be innacurate

%--------------TABLE GENERATION--------------
D=numel(LBs);
Ns = ones(1,D)*3;

dvars = {};
for dvar_idx = 1:numel(Ns)
    dvars{dvar_idx} = linspace(LBs(dvar_idx),UBs(dvar_idx),Ns(dvar_idx));
end

table = table_calculations(dvars,NaN,NaN,D,error_type,logspace_interpolation);

max_error = Inf;
idx = 0;
var = 0;
dvar_shifts = {};
table_shifts = {};
prev_dvar_expansion = 0;

total_time = 0;

while max_error > error_threshold
    %Determine errors associated with each input parameter individually
    for dvar_idx = 1:D
        dx = (dvars{dvar_idx}(2)-dvars{dvar_idx}(1))/2;
        dvar_shifts{dvar_idx} = linspace(dvars{dvar_idx}(1)+dx, dvars{dvar_idx}(end)-dx, numel(dvars{dvar_idx})-1); %#ok<*SAGROW>
        dvar_input = dvars;
        dvar_input{dvar_idx} = dvar_shifts{dvar_idx};
        if dvar_idx ~= prev_dvar_expansion && prev_dvar_expansion ~= 0
            dvar_input{prev_dvar_expansion} = dvar_input{prev_dvar_expansion}(2:2:end-1);
            [table_shifts_new{dvar_idx},errors_new{dvar_idx},max_errors(dvar_idx)] = table_calculations(dvar_input,dvars,table,D,error_type,logspace_interpolation);
            
            if D > 1
                dim_vector = linspace(1,D,D);
                dim_vector(prev_dvar_expansion) = 1;
                dim_vector(1) = prev_dvar_expansion;

                ptable_shift = permute(table_shifts{dvar_idx},dim_vector);
                ptable_shift_new = permute(table_shifts_new{dvar_idx},dim_vector);            

                table_shift_size = size(table_shifts{dvar_idx});
                table_shift_size(prev_dvar_expansion) = table_shift_size(prev_dvar_expansion)*2-1;

                table_shifts{dvar_idx} = permute(reshape([reshape([ptable_shift(1:end-1,:),ptable_shift_new(:,:)]',[size(ptable_shift(:,:),2),size(ptable_shift(1:end-1,:),1)+size(ptable_shift_new(:,:),1)])';ptable_shift(end,:)],table_shift_size(dim_vector)),dim_vector);

                perrors = permute(errors{dvar_idx},dim_vector);
                perrors_new = permute(errors_new{dvar_idx},dim_vector);  

                errors{dvar_idx} = permute(reshape([reshape([perrors(1:end-1,:),perrors_new(:,:)]',[size(perrors(:,:),2),size(perrors(1:end-1,:),1)+size(perrors_new(:,:),1)])';perrors(end,:)],table_shift_size(dim_vector)),dim_vector);
            
            else
                table_shift_size = size(table_shifts{1},1)*2-1;
                table_shifts{dvar_idx} = [reshape([table_shifts{dvar_idx}(1:end-1);table_shifts_new{dvar_idx}],[1,table_shift_size-1]),table_shifts{dvar_idx}(end)];
                errors{dvar_idx} = [reshape([errors{dvar_idx}(1:end-1);errors_new{dvar_idx}],[1,table_shift_size-1]),errors{dvar_idx}(end)];
                
            end

        else
            [table_shifts{dvar_idx},errors{dvar_idx},max_errors(dvar_idx)] = table_calculations(dvar_input,dvars,table,D,error_type,logspace_interpolation);
        end
    end
    

    %Approximate Errors
    [~,prev_dvar_expansion] = max(max_errors);
    if error_type ==0
        max_error = abs((1-prod(1+max_errors./100)))*100 %#ok<NOPTS>
    else
        max_error = sum(max_errors) %#ok<NOPTS>
    end
    
    
    %Add newly calculated values to the table
    if D > 1
        dim_vector = linspace(1,D,D);
        dim_vector(prev_dvar_expansion) = 1;
        dim_vector(1) = prev_dvar_expansion;
        
        ptable = permute(table,dim_vector);
        ptable_shift = permute(table_shifts{prev_dvar_expansion},dim_vector);   


        table_size = size(table);
        table_size(prev_dvar_expansion) = table_size(prev_dvar_expansion)*2-1;

        table = permute(reshape([reshape([ptable(1:end-1,:),ptable_shift(:,:)]',[size(ptable(:,:),2),size(ptable(1:end-1,:),1)+size(ptable_shift(:,:),1)])';ptable(end,:)],table_size(dim_vector)),dim_vector);
    
    else
        table_size = size(table,2)*2-1;
        table = [reshape([table(1:end-1);table_shifts{1}],[1,table_size-1]),table(end)];
    end

    dvars{prev_dvar_expansion} = linspace(LBs(prev_dvar_expansion),UBs(prev_dvar_expansion),table_size(prev_dvar_expansion));
    
    number_of_parameter_values = size(table)    
    
    idx = idx+1 %#ok<NOPTS>
end

%save table
dt = datestr(now, 'yyyymmdd_HH_MM_SS');
save(strcat(output_table_filepath,'\',dt,fn,'.mat'),'dvars','table')

%Function to calcualte linear interpolation errors associated with
%the specified numerical function
function [table_outputs,errors,max_error] = table_calculations(dvars,dvar_knowns,table_input,D,error_type,logspace_interpolation)
    
    size = zeros(1,D);
    dvar_vals_init = zeros(1,D);
    for dvar_idx = 1:D
        size(dvar_idx) = numel(dvars{dvar_idx});
        dvar_vals_init(dvar_idx) = dvars{dvar_idx}(1);
    end
    dvar_vals = dvar_vals_init;
    
    indicies = ones(1,D);
    
    table_outputs_1D_tot = [];
    errors_1D_tot = [];
    
    dvar_index_one = dvars{1};
    while true     
        %Below parallel process for loop (parfor) can be replaced by a
        %simple for loop if necessary.
        parfor idx = 1:size(1) 
            error_val = 0;
            temp = dvar_vals;
            temp(1) = dvar_index_one(idx);
            table_output = function_calculation(temp);           
            
            if ~isnan(table_input)
                dvar_vals_inputs = num2cell(temp);
                
                if logspace_interpolation == 1
                    approx = lininterpn(dvar_knowns{:},log(table_input),dvar_vals_inputs{:}); %#ok<*PFBNS>
                else
                    approx = lininterpn(dvar_knowns{:},table_input,dvar_vals_inputs{:});
                end
                if error_type == 0
                    if logspace_interpolation == 1
                        error_val = abs((table_output-exp(approx))/(abs(table_output)+abs(exp(approx))))*200;
                    else
                        error_val = abs((table_output-exp(approx))/(abs(table_output)+abs(approx)))*200;
                    end
                else
                    if logspace_interpolation == 1
                        error_val = abs(table_output-exp(approx));
                    else
                        error_val = abs(table_output-approx);
                    end
                end
            end        
            table_outputs_1D(idx) = table_output;
            errors_1D(idx) = error_val;
        end   
        dvar_vals(1) = dvars{1}(size(1));
        indicies(1) = size(1);
        
        table_outputs_1D_tot = [table_outputs_1D_tot,table_outputs_1D]; %#ok<*AGROW>
        errors_1D_tot = [errors_1D_tot,errors_1D];     
        
        [~,z] = find(indicies~=size,1);
        if isempty(z)
            break
        end
        
        indicies(z) = indicies(z) + 1;
        dvar_vals(z) = dvars{z}(indicies(z));
        indicies(1:z-1) = 1;
        dvar_vals(1:z-1) = dvar_vals_init(1:z-1);
    end
    
    table_outputs = reshape(table_outputs_1D_tot',size);
    errors = reshape(errors_1D_tot',size);
    
    max_error = max(errors(:));
end

%Function to solve specified numerical formula
function [calculation_output] = function_calculation(dvar_vals)
    phi = 4.8; %Work function (eV)
    [~,calculation_output] = FE_heat_exchange_numerical([dvar_vals,phi]);
end