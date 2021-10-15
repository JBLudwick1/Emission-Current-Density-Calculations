classdef NUMERICAL_FUNC
    methods ( Static = true )

        %----------Inputs----------
        %Fs:    Electric Field Magnitude [V/nm]
        %Ts:    Temperature     [K]
        %radii: Local Emitter Radius of Curvature [nm]
        %phis:  Work Function   [eV]
        %EFs:   Fermi Energy    [eV]

        %----------Outputs----------
        %Js:    Local Emission Current Density [A/cm^2]

        function Js = QMWI(Fs,Ts,radii,phis,EFs)

        global error_threshold;
        if isempty(error_threshold)
            error_threshold = 5;
        end

        % various constants needed in reduced units (i.e, normalized so we do not carry the large powers of ten)
        kB = 1.380649;         % normalized Boltzman Constant
        m0 = 9.1093837015;     % normalized Electron Effective Mass
        q = 1.602176634;       % normalized Electron Charge (magnitude)
        hbar = 1.054571817;    % normalized Planck's Constant /2pi
        h = 6.62607015;        % normalized Planck's Constant
        eps0 = 8.8541878128;   % normalized Vacuum Dielectric Constant
        PI = 3.14159265358979323846;

        A = (q*m0*kB^2)/(2*PI^2*hbar^3)*100;     %[A cm^-2 K^-2] Richardson-Laue-Dushman (RLD) Coeff
        B = 2*hbar/m0/1i*10^-3;

        %Count the length of all input arrays
        F_count = numel(Fs);
        T_count = numel(Ts);
        radii_count = numel(radii);
        phis_count = numel(phis);
        EFs_count = numel(EFs);

        %Initialize current density value matrix
        Js = zeros(T_count,F_count*radii_count*phis_count*EFs_count);

        %Used for determining energy range later
        kBT_max = kB*(3000)/q*1e-4;

        %Convert to a single parallel for loop and determine each parameter index from the
        %ndgrid variables
        lin = @(x) linspace(1,x,x);
        [F_indicies, radius_indicies, phi_indicies, EF_indicies] = ndgrid(lin(F_count), lin(radii_count), lin(phis_count), lin(EFs_count));

        parfor iteration = 1:numel(Js)/T_count

            % Obtain parameter indicies for iteration
            F_idx = F_indicies(iteration);
            radius_idx = radius_indicies(iteration);
            phi_idx = phi_indicies(iteration);
            EF_idx = EF_indicies(iteration);

            if any(isnan([Fs(F_idx), radii(radius_idx), phis(phi_idx), EFs(EF_idx)]))
                Js(:,iteration) = nan;
                continue;
            end

            EF = EFs(EF_idx);

            % Energy bounds to later calculate the transmission probability versus energy
            E_high = EF+30*kBT_max;
            E_low = 0;

            phi = phis(phi_idx);

            radius = radii(radius_idx);

            x0b = ((radius*((50*q)/PI + 4*EF*eps0*radius + 4*eps0*phi*radius)) / (eps0*(EF + phi)))^(1/2)/2 - radius;

            F = Fs(F_idx);

            % Determine 'Length': Size of the domain where the potential energy profile is defined
            % Set length to be equal to where E = Estar_factor*(-phi-EF). This length was tested
            % to be adequated between F = 0.1 to 10 V/nm and phi = 2 to 5 eV (radius = 5nm and 5m, EF = 5 eV)
            % -1e-6 term added to help with numerical instabilities when desired energy
            % is very near horizontal asymptote
            Estar_factor = 3;
            EpotFunc = @(x) EF + phi - F*radius*(x+x0b)./(radius + x + x0b) - (25*q*radius)./(2*eps0*pi*(x + x0b).*(2*radius + x + x0b));

            if phi+EF-F*radius <= Estar_factor*(-phi-EF)-1e-6        
                EpotFunc2 = @(x) EpotFunc(x) + Estar_factor*(phi+EF);
                Length = NUMERICAL_FUNC.bisection(EpotFunc2, 0, 200*radius);
                Emin = 0;
            else
                Length = 99*radius;
                Emin = phi+EF-F*radius;
            end

            if radius == Length %Fix error when radius = length
                Length = Length*1.01;
            end

            %% Discretize Potential Barrier into a series of step potentials
            % All the step functions are set to have a maximum allowed energy step so
            % that the step functions closely match the actual potential barrier
            % High slope barrier regions contain many step barriers over a short region of distance
            % Low slope barrier regions contain few step barriers over a lonng region of distance

            Epot = zeros(1,250);    % Barrier assumed to start at E = 0, x = 0
            xs = zeros(1,250);      % Initialize x array
            max_dE = 0.1;           % Max step potential height
            E_idx = 2;              % Since Epot is known at x=0, start at index = 2
            dx = 0.05;              % Inital distance step guess
            x = 0;

            while x < Length
                x = xs(E_idx-1) + dx;               % Set x location guess
                Epotval = EpotFunc(x);

                % Check if step barrier height is lower than max threshold (max_dE)
                if ~(abs(Epot(E_idx-1) - Epotval) <= max_dE || dx < 1e-6)
                    % Step barrier height higher than maximum threshold (max_dE)
                    % x location guess must be smaller
                    % Reduce x step by a factor of 90%
                    % Reducing this value too quickly may result in an excessive number of step barriers
                    % Reducing too little increasing the computational toll of this while loop
                    dx = dx*0.9;
                    continue;
                end

                dx = abs(x - xs(E_idx-1))*1.1;      % next init x step guess
                xs(E_idx) = x;                      % Add step function position to x array
                Epot(E_idx) = Epotval;              % Add step function energy to Epot array
                E_idx = E_idx + 1;                  % Incrememnt index
            end

            xs(E_idx:end) = [];
            Epot(E_idx:end) = [];

            Epot = (Epot(2:end) + Epot(1:end-1))/2;
            xs = xs(1:end-1);
            dxs = xs(2:end) - xs(1:end-1);                

            %% Determine sufficient energy step resolution for integration by trapezoidal
            % method by doubling energy resolution and monitoring percent error change
            Esteps = 100;            
            Es = linspace(E_low,E_high,Esteps);
            Es_tot = Es;

            N = zeros(numel(Es), T_count);
            J = zeros(1,T_count);

            error = 100;
            R_old = [];
            J_old = NaN;
            Rflag = true;
            while error > error_threshold/2 || any(isnan(error))
                %% Vectorized
                E = Es;

                % QMWI APPROACH------------------------------------------
                k_inc = sqrt(2*m0*q*E)/hbar;    

                % iterating over domain to find Z(0+) for a specific incident energy      
                scalar = 1 / hbar * sqrt(2 * q * m0);
                ZR = 2e6 * sqrt(2/m0*q*complex(E-Epot(end)));
                ZR = arrayfun( @(x,y) NUMERICAL_FUNC.CalcImpedance2(x,y,dxs,Epot,scalar), (ZR), (E));

                % Final Z value is Z(0+)
                Z0plus = ZR;
                Z0 = 2*k_inc*hbar/m0*10^6;

                % Calculate transmission probability at specific incident energy   
                % Solve the transmission probability with the reflection coefficient.
                % This works well at T near 1, but not as T approaches zero
                rho = (Z0 - Z0plus)./(Z0 + Z0plus);
                Trans1 = (1-abs(rho).^2);

                % Solve the transmission probability with the transmission coefficient.
                % This works well at low T, but not as T approaches unity
                t = 2*Z0./(Z0+Z0plus);
                Trans2 = abs(t).^2.*real(Z0plus./Z0);
                if(E(1) == 0), Trans2(1) = 0; end

                test = find(Trans1<=1e-4,1,'last');
                if isempty(test)
                    Trans = Trans1;
                else
                    Trans = [Trans2(1:test), Trans1(test+1:end)];
                end

                % Determine supply function values
                [E2, Ts2] = ndgrid(E,Ts);

                reg1 = exp(-q*(E2 - EF)./(kB*Ts2)*10^4) < 1e-5;                   %Use approximate solutions to avoid numerical precision errors
                reg2 = exp(-q*(E2 - EF)./(kB*Ts2)*10^4) > 1e5;                    %Use approximate solutions to avoid numerical precision errors
                reg3 = ~(reg1 | reg2);

                N(reg1) = exp(-q*(E2(reg1) - EF)./(kB*Ts2(reg1))*10^4);           %Use approximate solutions to avoid numerical precision errors
                N(reg2) = -q*(E2(reg2) - EF)./(kB*Ts2(reg2))*10^4;                %Use approximate solutions to avoid numerical precision errors
                N(reg3) = log(exp(-q*(E2(reg3) - EF)./(kB*Ts2(reg3))*10^4) + 1);   %Direct formula of N in a region where numerical precision is fine

                Trans(E < Emin | E == 0) = 0;

                [Trans2, ~] = ndgrid(Trans,Ts);
                R = N .* Trans2; %Emission Current Density per normal energy            

                %% Expand energy resolution
                if numel(R_old) > 0
                    R_new = zeros(numel(Es_tot),numel(Ts));
                    R_new(1:2:end,:) = R_old;
                    R_new(2:2:end-1,:) = R;
                    R = R_new;
                end                   

                % integrate and calculate J
                den = trapz(Es_tot,R);
                J = A*Ts/kB*q.*den*10^4; % Calculate Emission Current Density

                error = max(abs((J-J_old)./J))*100; % Compare emission current density percent change when energy resolution is doubled

                if Rflag
                    Rflag = false;
                    percent = 0.1/100;

                    % Find left R value that is less than 0.1 percent of the left most peak
                    [M,I] = max(R(:,1));
                    LeftR = find(R(1:I,1) < M*percent,1,'last');
                    if isempty(LeftR), LeftR = 1; end

                    % Change the lower E bound
                    LeftR = max(1,LeftR-1);
                    if E(LeftR) > E_low
                    E_low = E(LeftR);
                    end

                    % Find right R value that is less than 0.1 percent of the right most peak
                    [M,I] = max(R(:,end));
                    RightR = find(R(I:end,end) < M*percent,1,'first') + I - 1;
                    if isempty(RightR), RightR = size(R,1); end

                    % Change the higher E bound
                    RightR = min(size(R,1),RightR+1);
                    if E(RightR) < E_high
                    E_high = E(RightR);
                    end

                    % Focus the steps in the area where R is
                    Esteps = RightR - LeftR + 1;
                    % Esteps = min(Esteps,50);
                    R = R(LeftR:RightR,:);
                end



                % Create array of new energy values to increase resolution
                Esteps = Esteps * 2 - 1;
                Es_tot = linspace(E_low,E_high,Esteps);
                Es = Es_tot(2:2:end);

                J_old = J; % Initialize for next iteration
                R_old = R; % Initialize for next iteration
                R = zeros(numel(Es), T_count);
                N = zeros(numel(Es), T_count);
            end

            % Use final value of emission current density with
            % sufficient energy resolution
            Js(:,iteration) = J;
        end

        % Ensure table output form is what lookup table code expects
        Js = reshape(Js,[T_count,F_count,radii_count,phis_count,EFs_count]); 

        % Reorder to have temperature as dimension 2, and field as dimension 1 
        Js = permute(Js,[2,1,3,4,5]);
        end

        function ZR0 = CalcImpedance2(ZR0, E0, dxs, Epot, scalar)
            % iterating over domain to find Z(0+) for a specific incident energy
            B = -2.315352721010859e-04i;
            for Epidx=numel(dxs):-1:1
                sqrtval = sqrt(abs(E0-Epot(Epidx)));
                k = scalar * sqrtval * dxs(Epidx);
                Bk = -scalar * sqrtval * B * 1e9;

                if E0 > Epot(Epidx)
                    ZR0 = Bk * (ZR0 * cos(k)  - Bk * sin(k))  / (ZR0 * sin(k)  + Bk * cos(k));
                else
                    ZR0 = Bk * (ZR0 * cosh(k) + Bk * sinh(k)) / (ZR0 * sinh(k) + Bk * cosh(k)); 
                end
            end
        end 

        function c = bisection(func, a, b)
            c = a;
            EPSILON = 1e-10;
            while (b - a) >= EPSILON
                c = (a + b) / 2;                    % Find middle point

                if func(c) == 0; break; end         % Check if middle point is root

                if sign(func(c) * func(a)) == 1     % Decide the side to repeat the steps
                    a = c;
                else
                    b = c;
                end
            end
        end
    end
end