classdef linear_interp
    methods ( Static = true )
		% linear interpolation - input 1-dimensional arrays X1, X2, ... , Xn, a n-dimensional array V, and x1, x2, ..., xn query values
		% assumes Xi's are in increasing order. Can query multiple points, but x1,
		% x2, ..., xn must be the same size
		%
		% Differences from matlab built-in :
		%       much, much faster
		%       if coordinate is exactly on the spot, doesn't look at neighbors.  e.g. interpolate([blah, blah2], [0, NaN], blah) returns 0 instead of NaN
		%       extends values off the ends instead of giving NaN
		%

		%This function relies on X and Y to be a strictly increasing vectors
		function v = interpn(varargin)         
            if (mod(numel(varargin), 2) ~= 1), error('Invalid number of inputs.  Should be odd'); end
            n = (length(varargin) - 1) / 2;
            varargin{n+1} = squeeze(varargin{n+1});
            N = size(varargin{n+1});
            N(N==1) = [];

            if any( cellfun(@numel, varargin(1:n)) ~= N)
                error('length of parameters do not match size of table');
            elseif ~all( cellfun(@numel, varargin(n+2:end)) == numel(varargin{n+2}) )
                error('length of query values are not equal');
            end

            % Convert all to columns
            varargin([1:n,n+2:end]) = cellfun(@(c) reshape(c,[],1), varargin([1:n,n+2:end]),'un',false);
            if isrow(varargin{n+1}); varargin{n+1} = varargin{n+1}.'; end

            pindex = cell(n, 1);
            oindex = cell(n, 1);
            slopes = cell(n, 1);

            for i = 1:1:n
                x = varargin{n+1+i};
                X = varargin{i}.';
                [pindex{i}, oindex{i}, slopes{i}] = linear_interp.FSlope(x, X, numel(X));
            end

            V = varargin{n+1};
            v = zeros(size(varargin{end}));
            
            for bin = 1 : 1 : 2^n
                indexgetter = bin;
                multiplier = ones(numel(pindex{i}), 1);
                indices = cell(n, 1);

                for i = 1:1:n
                    index = mod(indexgetter, 2);
                    indexgetter = (indexgetter - index)/2;
                    if index == 0 
                        indices{i} = pindex{i};
                        multiplier = multiplier .* (1 - slopes{i});
                    else
                        indices{i} = oindex{i};
                        multiplier = multiplier .* slopes{i};
                    end
                end

                lidx = sub2ind(size(V), indices{:});
                v = v + V(lidx) .* multiplier;                
            end

        end

        function [pindex, oindex, slope] = FSlope(a, A, ALen) 
            %Vector version of Fslope
            %Requires that x query values are column ie [1; 2; 3] and y query values
            %are row [1 2 3]
            
            pindex = zeros(size(a));
            oindex = zeros(size(a));
            slope = zeros(size(a));

            %Where a is less than A(1)
            less = a <= A(1);
            pindex(less) = 1;
            oindex(less) = 1;

            %Where a is greater than the end of A
            greater = a >= A(end);
            pindex(greater) = ALen;
            oindex(greater) = ALen;

            %Where it is in between the bounds
            %May need to add an epsilon value for checking,
            %but I think matlab handles it fairly well
            neither = ~(less | greater);
            
            temp = a(neither);
            
            %Check if query value directly at table value
            identicalcheck  = false(numel(a),1);
            
            if ~isempty(temp)
                [~,~,B]= histcounts(temp,A);
                oindex(neither) = B;
                pindex(neither) = B + 1;
                
                %Check if query value directly at table value
                ErrorMax = 1e-7;
                identicalcheck (neither) = abs((temp-A(B).')./temp)<ErrorMax;

                pindex(identicalcheck ) = B(identicalcheck (neither));
            end
            
            %Not directly on table value and inbetween the bounds
            neither = (neither & ~identicalcheck );

            %Calculate slope
            Xp = zeros(size(a));
            Xi = zeros(size(a));

            Xp(neither) = A(pindex(neither));
            Xi(neither) = A(oindex(neither));

            slope(neither) = (a(neither) - Xp(neither)) ./ (Xi(neither) - Xp(neither));
        end
	end
end