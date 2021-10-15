classdef LUT_gen
    methods ( Static = true )
        function [] = lookup_table_generator(LBs_all,UBs_all,special_funcs,special_funcs_inv,logspace_interpolation,output_full_path,num_fun_str)
        global error_threshold;
        if isempty(error_threshold)
            error_threshold = 5;
        end
        
        fprintf('Output path: %s\n', output_full_path);
        fprintf('Lower bounds: %s\n', sprintf('%8.2f ', LBs_all));
        fprintf('Upper bounds: %s\n', sprintf('%8.2f ', UBs_all));

        %--------------TABLE GENERATION--------------
        %Find constants where lower and upper bounds are the same
        constants = UBs_all(UBs_all==LBs_all);
        
        %Convert specified parameters using special functions
        for idx = 1:numel(UBs_all)
            if ~isnumeric(special_funcs{idx})
                UBs_all(idx) = special_funcs{idx}(UBs_all(idx));
                LBs_all(idx) = special_funcs{idx}(LBs_all(idx));
                fprintf('Special function for parameter %s: %s\n', num2str(idx),func2str(special_funcs{idx}));
            end
        end
        
        %Remove constants from sweep array
        dvar_indicies = find(UBs_all~=LBs_all);
        UBs = UBs_all(dvar_indicies);
        LBs = LBs_all(dvar_indicies);
        
        D = numel(UBs); %Determine lookup table dimensionality
        Ns = ones(1,D)*3; %Inital set size of each input parameter to 3 values each

        fprintf('Error threshold: %5d\n', error_threshold);
        fprintf('Logspace interp?: %s\n', num2str(logspace_interpolation));
        fprintf('Starting dim: %s\n\n', sprintf('%8d ', Ns));

        % dvars contains all the input parameter values
        dvars = cell(1,D);
        dvarsExpMid = cell(1,D);
        dvarsExpInterMid = cell(1,D);

        for dvar_idx = 1:D
            dvars{dvar_idx} = linspace(LBs(dvar_idx),UBs(dvar_idx),Ns(dvar_idx));
            dvarsExpMid{dvar_idx} = conv(double(true(size(dvars{dvar_idx}))), [0.5 0.5], 'valid');
            dvarsExpInterMid{dvar_idx} = conv(double(true(size(dvars{dvar_idx}))), [0.5 0.5], 'valid');
        end

        dvarsExpMidNew = cell(1,D);
        dvarsExpInterMidNew = cell(1,D);
        dvar_table_shift = cell(1,D);
        dvar_shifts_previous = cell(1,D);
        checkidx = cell(1,D);

        % Initial table calculation performed. NaN inputs indicate no error
        % measurements are being conducted here
        table = LUT_gen.table_calculations(dvars,NaN,NaN,logspace_interpolation,num_fun_str,constants,dvar_indicies,special_funcs,special_funcs_inv);

        % Initialization of variables
        max_error = Inf;
        idx = 0;
        prev_dvar_expansion = 0;

        dvar_shifts = cell(1,D);
        table_shifts = cell(1,D);
        table_shifts_new = cell(1,D);
        errors = cell(1,D);
        errors_new = cell(1,D);
        max_errors = zeros(1,D);

        % Start total time timer
        totaltime = tic;

        % While loop to repeat growing lookup table until the maximum predicted
        % error drops below a desired error threshold
        while max_error > error_threshold
            timestart = tic; % Start iteration timer

            for dvar_idx = 1:D % Loop through each parameter
%                 titlewb = sprintf('Please wait... Calculating for variable #%d',dvar_idx);

                % Expand it how we calculated earlier by removing extra terms
                if dvar_idx ~= prev_dvar_expansion && prev_dvar_expansion ~= 0
                    checkidx{dvar_idx} = dvarsExpMid{dvar_idx};
                else
                    checkidx{dvar_idx} = dvarsExpInterMid{dvar_idx};
                end

                % If this dvar doesn't have the worst error and we have already
                % calculated some values for it, expand upon those already calculated
                % values
                if dvar_idx ~= prev_dvar_expansion && prev_dvar_expansion ~= 0
                    % Grab the expanded values of the dvar that expanded last, and continue
                    % in that direction as we've calculated all the other terms
                    dvar_shifts{dvar_idx} = dvar_shifts_previous{dvar_idx};
                    dvar_input = dvars;

                    if numel(dvar_shifts_previous{dvar_idx}) ~= sum(checkidx{dvar_idx})
                        dvar_shifts_previous{dvar_idx}(~checkidx{dvar_idx}) = [];
                    end

                    dvar_input{dvar_idx} = dvar_shifts_previous{dvar_idx};
                    dvar_input{prev_dvar_expansion} = dvar_shifts_previous{prev_dvar_expansion};

                    if isempty(dvar_input{dvar_idx})
                       dvar_table_shift{dvar_idx} = dvar_input{dvar_idx};
                       continue 
                    end

                    [table_shifts_new{dvar_idx},errors_new{dvar_idx},max_errors(dvar_idx)] = LUT_gen.table_calculations(dvar_input,dvars,table,logspace_interpolation,num_fun_str,constants,dvar_indicies,special_funcs,special_funcs_inv);

                    if D > 1
                        % Remove any columns that we thought we were going to use,
                        % and now don't, before permuting dimensions
                        if size(table_shifts{dvar_idx},dvar_idx) ~= size(table_shifts_new{dvar_idx},dvar_idx)
                            % Store what dvars were used to make this
                            % table, and then compare that to the dvars I want to
                            % use now, and delete the extra
                            dim_vector = linspace(1,D,D);
                            dim_vector(dvar_idx) = 1;
                            dim_vector(1) = dvar_idx;

                            check2 = sum(repmat(dvar_table_shift{dvar_idx},[numel(dvar_input{dvar_idx}),1]) == repmat(dvar_input{dvar_idx}',[1,numel(dvar_table_shift{dvar_idx})]),1);

                            ptable = permute(table_shifts{dvar_idx},dim_vector);
                            perrors = permute(errors{dvar_idx},dim_vector);

                            ptable_size = size(ptable); 
                            ptable_size(1) = ptable_size(1) - sum(~check2);
                            perrors_size = size(perrors); 
                            perrors_size(1) = perrors_size(1) - sum(~check2);

                            ptable(~check2,:) = []; 
                            perrors(~check2,:) = []; 

                            ptable = reshape(ptable,ptable_size);
                            perrors = reshape(perrors,perrors_size);

                            table_shifts{dvar_idx} = permute(ptable,dim_vector);
                            errors{dvar_idx} = permute(perrors,dim_vector);
                        end

                        % Reorder the dimensions so the arrays that we want to
                        % concatenate are in the first dimension
                        dim_vector = linspace(1,D,D);
                        dim_vector(prev_dvar_expansion) = 1;
                        dim_vector(1) = prev_dvar_expansion;

                        % Reorder the dimensions of the old and new shifted tables
                        % and errors matrices so they can be concatenated, in order of dim_vector
                        ptable_shift = permute(table_shifts{dvar_idx},dim_vector);
                        ptable_shift_new = permute(table_shifts_new{dvar_idx},dim_vector);
                        perrors = permute(errors{dvar_idx},dim_vector);
                        perrors_new = permute(errors_new{dvar_idx},dim_vector);

                        % Concatenate the dvars and then sort them so we know how to sort
                        % the newly concatenated but out of order table
                        concat_dvars = [dvar_previous{prev_dvar_expansion}, dvar_shifts_previous{prev_dvar_expansion}];
                        [~,I] = sort(concat_dvars); 

                        % Concatenate the tables
                        table_shifts{dvar_idx} = cat(1, ptable_shift, ptable_shift_new);
                        table_shift_size = size(table_shifts{dvar_idx});

                        % Reorder the newly made table in order of I
                        % This assumes that dvars are in a linearly increasing order normally.
                        % As well as undo the dimension shifting done earlier
                        table_shifts{dvar_idx} = table_shifts{dvar_idx}(I,:);
                        table_shifts{dvar_idx} = reshape(table_shifts{dvar_idx},table_shift_size);
                        table_shifts{dvar_idx} = permute(table_shifts{dvar_idx},dim_vector);

                        % Concatenate the errors
                        errors{dvar_idx} = cat(1, perrors, perrors_new);
                        error_shift_size = size(errors{dvar_idx});

                        % Reorder the newly made errors in order of the dvars and
                        % undo the dimension shifting done earlier
                        errors{dvar_idx} = errors{dvar_idx}(I,:);
                        errors{dvar_idx} = reshape(errors{dvar_idx},error_shift_size);
                        errors{dvar_idx} = permute(errors{dvar_idx},dim_vector);
                    else
                        table_shift_size = size(table_shifts{1},1)*2-1;
                        table_shifts{dvar_idx} = [reshape([table_shifts{dvar_idx}(1:end-1);table_shifts_new{dvar_idx}],[1,table_shift_size-1]),table_shifts{dvar_idx}(end)];
                        errors{dvar_idx} = [reshape([errors{dvar_idx}(1:end-1);errors_new{dvar_idx}],[1,table_shift_size-1]),errors{dvar_idx}(end)];
                    end
                else
                    % Create parameter array where each value is exactly halfway
                    % between each adjacent dvar points
                    dvar_shifts{dvar_idx} = conv(dvars{dvar_idx}, [0.5 0.5], 'valid');

                    % Grab all of the other dvars, and then calculate a sub table
                    % used to shift our final table
                    dvar_input = dvars; 
                    dvar_input{dvar_idx} = dvar_shifts{dvar_idx};
                    dvar_input{dvar_idx}(~checkidx{dvar_idx}) = []; 

                    if isempty(dvar_input{dvar_idx})
        %                continue 
                    end

                    [table_shifts{dvar_idx},errors{dvar_idx},max_errors(dvar_idx)] = LUT_gen.table_calculations(dvar_input,dvars,table,logspace_interpolation,num_fun_str,constants,dvar_indicies,special_funcs,special_funcs_inv);
                end

                % Save the values actually used for shifting
                dvar_table_shift{dvar_idx} = dvar_input{dvar_idx};
            end

            if (prev_dvar_expansion == 0); errors_new = errors; end

            % Find the maximum error for each dvar point used, and use that to reduce
            % the number of computations, by making a vector mask
            for i = 1 : D
                %Look at all dimensions except the current dim
                j = 1 : D; j(j==i) = []; 

                if isempty(dvar_table_shift{i})
                    continue;
                end

                % Expand the worst error dvar 
                if (i == prev_dvar_expansion)
                    [dvarsExpInterMidNew{i}, ~, dvarsExpMidNew{i}, ~] = LUT_gen.DvarExpansionMask(errors{i}, j, dvarsExpInterMid{i});
                else
                    [dvarsExpInterMidNew{i}, ~, dvarsExpMidNew{i}, ~] = LUT_gen.DvarExpansionMask(errors{i}, j, dvarsExpMid{i});
                end
            end

            % Determine parameter corresponding to the highest error
            [~,prev_dvar_expansion] = max(max_errors);              

            % Max value out of estimated 2D errors and 1D errors
            % The 2D should mostly always be worse
            if sum(isnan(max_error)) > 0
                max_error = NaN;
            else
                max_error = max(max_errors);   
            end

            if max_error > error_threshold
                % Reorder the dimensions so the arrays that what we want to
                % concatenate are in the first dimension
                dim_vector = linspace(1,D,D);
                dim_vector(prev_dvar_expansion) = 1;
                dim_vector(1) = prev_dvar_expansion;

                % Reorder the old table and the newly made table values
                % so they can be concatenated, in order of dim_vector
                if D > 1
                    ptable = permute(table,dim_vector);
                    ptable_shift = permute(table_shifts{prev_dvar_expansion},dim_vector);
                else
                    ptable = table;
                    ptable_shift = table_shifts{prev_dvar_expansion};
                end

                % Concatenate the tables
                table = cat(1, ptable, ptable_shift);
                tablesize = size(table);

                % Concatenate the dvars so we know how to sort the newly concatenated but out of order table
                % Sort the dvars that made this matrix in increasing order
                prev_dvar_shift = dvar_shifts{prev_dvar_expansion};
                prev_dvar_shift(~checkidx{prev_dvar_expansion}) = [];

                concat_prev_dvar = [dvars{prev_dvar_expansion}, prev_dvar_shift];
                [concat_prev_dvar,I] = sort(concat_prev_dvar); 

                % Save these values
                dvar_shifts_previous = dvar_shifts;
                dvar_shifts_previous{prev_dvar_expansion} = prev_dvar_shift;
                dvar_previous = dvars;

                % Reorder the newly made table in order of I
                % This assumes that dvars are in a linearly increasing order normally.
                % As well as undo the dimension shifting done earlier
                table = table(I,:);
                table = reshape(table,tablesize());
                if D > 1
                    table = permute(table,dim_vector);
                end
                 % Save the dvars
                dvars{prev_dvar_expansion} = concat_prev_dvar;
            end

            %Update how to expand values
            dvarsExpMid = dvarsExpMidNew;
            dvarsExpInterMid = dvarsExpInterMidNew;

            %Check current size of lookup table
            sizes = size(table);   

            %Increment index
            idx = idx+1;

            %Print Progress
            index_time = toc(timestart) / 60;
            total_time = toc(totaltime) / 60;

            fprintf('Iteration: %d\n',idx)
            fprintf('Approximated Maximum Error: %.2f%%\n',max_error)
            for pidx = 1:D
                fprintf('Parameter: %d, size: %d, max error: %.2f%%\n',pidx,sizes(pidx),max_errors(pidx))
            end
            fprintf('This Iteration Run Time (min): %.2f\n',index_time)
            fprintf('Total Run Time (min): %.2f\n\n',total_time)  
            
            idx_offset = 0;
            for special_idx = 1:numel(special_funcs_inv)
                if UBs_all(special_idx) == LBs_all(special_idx)
                    idx_offset = idx_offset+1;
                elseif ~isnumeric(special_funcs_inv{special_idx-idx_offset})
                    dvars_linspace{special_idx-idx_offset} = special_funcs_inv{special_idx-idx_offset}(dvars{special_idx-idx_offset}); 
                else
                    dvars_linspace{special_idx-idx_offset} = dvars{special_idx-idx_offset};
                end                
            end

            %save table (saves each iteration in case of runtime failure)            
            save(output_full_path,'error_threshold','max_errors','dvars_linspace','table','total_time','logspace_interpolation', '-v7.3')

        end
        end
   
        function [table_outputs,errors,max_error] = table_calculations(dvars,dvar_knowns,table_input,logspace_interpolation,num_fun_str,constants,dvar_indicies,special_funcs,special_funcs_inv)

            table_outputs = LUT_gen.function_calculation(dvars,num_fun_str,constants,dvar_indicies,special_funcs,special_funcs_inv);

            %Get values directly from known table values (where there is no error)
            numerical = table_outputs(:);

            if logspace_interpolation == 1
                table_outputs = log(table_outputs);
            end

            Ns = cellfun('prodofsize', dvars);
            errors = zeros([Ns,1]);

            if numel(table_input) > 1        
                % Get every combination of dvars
                [dvar_int{1:numel(dvars)}]=ndgrid(dvars{:});
                dvar_int = cellfun(@(c) reshape(c,[],1), dvar_int,'UniformOutput',false);

                %Approximate values directly between 2 known table values (where error is assumed to be worst) using linear interpolation 
                approx = linear_interp.interpn(dvar_knowns{:}, table_input, dvar_int{:});

                %Convert out of logspace for error calculations
                if logspace_interpolation == 1
                    approx = exp(approx);
                end

                %Calculate error values directly between 2 known table values (where error is assumed to be worst)
                errors = abs((numerical-approx)./(numerical+approx))*200;

                %Avoid divide by zero errors
                errors(numerical+approx == 0) = 0;

                %Catch any potential numerical errors here and assume a zero error
                errors(~isfinite(errors)) = 0;
            end

            %Find the max error values associated with each parameter
            max_error = max(errors(:));

            if isempty(max_error)
                max_error = 0;
            end

            if numel(dvars) > 1
                errors = reshape(errors, Ns);
            end
        end
        
        function [dvarsExpInterMid, dvarsExpInter, dvarsExpMid, dvarsExp] = DvarExpansionMask(Matrix, Dim, MidVectorUsed)
            %----------Inputs----------
            % Matrix:           the error matrix made by all the middle values of the current dvars
            % Dim:              the dimensions that MidVectorUsed does not expand
            % MidVectorUsed:    the midpoints used to generate the error Matrix

            %----------Outputs----------
            % We need to output all of the possible values that the dvar may use,
            % whether it is the worst error dvar or not    

            global error_threshold;
            if isempty(error_threshold)
                error_threshold = 5;
            end

            % What are the locations of the midpoints that were used
            MidIndices = false(1,numel(MidVectorUsed)*2+1);
            MidIndices(2:2:numel(MidVectorUsed)*2) = MidVectorUsed;

            % Used to delete some vals later
            UnusedMidElements = false(1,numel(MidVectorUsed)*2+1);
            UnusedMidElements(2:2:numel(MidVectorUsed)*2) = logical(MidVectorUsed==0);

            Midpoints = true(1,sum(MidVectorUsed));

            % Find the max error due to the full and mid points
            % If max errors fall below min threshold, set to zero
            maxErrors = Matrix;
            for didx = numel(Dim):-1:1
                maxErrors = max(maxErrors,[],Dim(didx));
            end
            Midpoints(maxErrors < error_threshold/2) = false;

            % Set indices to point values according to threshold
            MidIndices(MidIndices~=0) = Midpoints;
            MidVectorUsed2 = MidIndices(2:2:end);

            FullVectorUsed = false(1,numel(MidVectorUsed)+1);
            for idx = 1 : numel(MidVectorUsed2)
               if MidVectorUsed2(idx)
                   FullVectorUsed(idx:idx+1) = true;
                   continue;
               end
               FullVectorUsed(idx+1) = false;
            end

            % What are the locations of the fullpoints that were used
            FullIndices = zeros(1,numel(FullVectorUsed)*2-1);
            FullIndices(1:2:numel(FullVectorUsed)*2) = FullVectorUsed;

            % Inverleave the vector in case 
            Interleave = MidIndices + FullIndices;
            Interleave(UnusedMidElements) = [];
            FullIndices(2:2:end) = [];
            MidIndices(1:2:end) = [];

            % Figure out what mid points to use of the new interleaved values
            InterleaveMid = logical(floor(conv(Interleave, [0.5 0.5], 'valid')));

            dvarsExp = logical(FullIndices);
            dvarsExpMid = logical(MidIndices);
            dvarsExpInter = logical(Interleave);
            dvarsExpInterMid = logical(InterleaveMid);
        end

        function [calculation_output] = function_calculation(dvars,num_fun_str,constants,dvar_indicies,special_funcs,special_funcs_inv)

            if any(cellfun(@(C) any(isnan(C(:))), dvars)) || ...
                any(cellfun(@(C) any(isempty(C(:))), dvars))
                calculation_output = [];
                return;
            end
            
            %Arrange inputs to consist of sweep parameters and inputs
            N = numel(dvars)+numel(constants);
            inputs = cell(1,N);
            const_idx = 1;
            dvar_idx = 1;
            for idx = 1:N
                if sum(dvar_indicies == idx) > 0
                    inputs{idx} = dvars{dvar_idx};
                    dvar_idx = dvar_idx+1;
                else                    
                    inputs{idx} = constants(const_idx);
                    const_idx = const_idx+1;
                end
                if ~isnumeric(special_funcs{idx}); inputs{idx} = special_funcs_inv{idx}(inputs{idx}); end
            end
            
            %Call numerical function
            num_fun = str2func(num_fun_str);
            [calculation_output] = num_fun(inputs{:});
            
            %Reorder matrix
            const_indicies = linspace(1,N,N);
            const_indicies(dvar_indicies) = [];
            dim_order = [dvar_indicies,const_indicies];
            calculation_output = permute(calculation_output,dim_order);
        end
    end
end