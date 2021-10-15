function output = MATLAB_interpolate(F_vals,T_vals,R_vals,phi_vals)

    %WARNING: this code does not catch when query values are outside the
    %LUT bounds. If this happens, the query value becomes the closest valid
    %value in the LUT.

    %Load Lookup table. This global variable ensures the table is only
    %loaded one time
    global LUT
    if isempty(LUT)
        LUT = load('LUT_4D.mat','dvars_linspace','table');
    end
    %LUT is in logsapce so 'exp' converts to current density in A/cm^2
    %1e4 scalar converts to A/m^2
    output = exp(linear_interp.interpn(LUT.dvars_linspace{:},LUT.table,F_vals,T_vals,R_vals,phi_vals)')*1e4;
    
    %if any electric field values fall below the minimum LUT range, assume
    %zero current
    output(F_vals<min(LUT.dvars_linspace{1})) = 0;
    
end