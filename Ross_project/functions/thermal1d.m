classdef thermal1d
    % thermalSolution1d Implements a numerical 1D heat transfer model for an ice sheet
    % Author: Donglai Yang, dyang379@gatech.edu
    % Sep 2024
    
    properties
        % Public properties that can be accessed directly
        thickness   % Thickness of the ice sheet
        temperature % Temperature profile of the ice sheet
        dsdx % surface slope
        z % vertical axis
        smb % surface mass balance (m/a); positive is accumulating
        bm % basal melt rate (m/a); positive is melting
        Ts % surface temperature (K)
        G % geothermal heat flux (W/m^2)
        u % along flow velocity (m/a)
        dTdx % along flow temperature gradient (K/m)
        w_flag % flag for vertical velocity formulation: "linear" or "DJ" (Dansgaard-Johnsen) or "Nye"
        solver % solver (either collocation method or finite difference)
        sol % Temperature solution (K)
    end
    
    properties (Access = private)
        % Private properties, only accessible within the class
        dz          % Spatial step size
    end
    
    properties (Constant) % accessing these uses <classname>.<constant_name>, not obj.<constant_name>
        % as it is associated with the entire class, not an particular
        % instance
        % Constants that all methods can access
        k = 2.22;  % W/(m·K) thermal conductivity for ice
        cp = 2050;       % J/(kg·K) heat capacity for ice
        rho_i = 917;     % kg/m^3 ice dnesity
        n = 3;  % unitless, Glen's flow law constant
        g = 9.8; % OH GRAVITY
        Tm = 273.15; % K, melting point
        L = 1e3; % characteristic glacier length (shouldn't matter)
    end
    
    methods % dynamic methods, requiring obj.
        function obj = thermal1d(thickness, slope, dz, smb, bm, Ts, G, u, dTdx, w_flag, solver)
            % Constructor: Initialize the thermal solution
            obj.thickness = thickness;
            obj.dsdx = slope;
            obj.dz = dz;
            obj.z = thickness:-dz:0; % surface to base
            obj.smb = smb;
            obj.bm = bm;
            obj.Ts = Ts; 
            obj.G = G;
            obj.u = u;
            obj.dTdx = dTdx;
            obj.w_flag = w_flag;
            obj.solver = solver;
        end

        function dydx = thermalLinearSystem(obj, z, y)
            % represent 1D heat transfer as a system of 1st order ODE
            % z is the independent variable (depth)
            % y(1) is T (temperature); y(2) is temperature gradient (dT/dz)

            T_assume = 260; % assumed temperature for calculating prefactor A

            switch obj.w_flag
                case "linear"
                    if obj.smb == 0; disp('SMB given = 0, hence w=0!');end
                    w = obj.calcLinearVerticalVel(z);
                case "DJ"
                    w = obj.calcDJVerticalVel(z);
                case "Nye"
                    w = obj.calcNyeVerticalVel(z);
                case "constant"
                    w = -thermal1d.yrToSec(obj.smb + obj.bm);
                otherwise
                    error("Vertical velocity request unknown!")
            end
            Q = obj.calcViscousDissipationFromT(z, T_assume);
            
            dydx = zeros(2,1);
            dydx(1) = y(2);
            dydx(2) = y(2).*w./thermal1d.getAlpha()...
                      - Q./(thermal1d.rho_i*thermal1d.cp*thermal1d.getAlpha());
        end

        function dydx = thermalLinearSystemHAC(obj, z, y)
            % represent 1D heat transfer as a system of 1st order ODE
            % z is the independent variable (depth)
            % y(1) is T (temperature); y(2) is temperature gradient (dT/dz)
            % HAC stands for with horizontal advective cooling; here we
            % also ignore viscous disspation since it is used for ice shelf
            % (no basal drag and hence it is a plug flow)

            switch obj.w_flag
                case "linear"
                    if obj.smb == 0; disp('SMB given = 0, hence w=0!');end
                    w = obj.calcLinearVerticalVel(z);
                case "DJ"
                    w = obj.calcDJVerticalVel(z);
                case "Nye"
                    w = obj.calcNyeVerticalVel(z);
                case "constant"
                    w = -thermal1d.yrToSec(obj.smb + obj.bm);
                otherwise
                    error("Vertical velocity request unknown!")
            end
            
            dydx = zeros(2,1);
            dydx(1) = y(2);
            dydx(2) = 1./thermal1d.getAlpha().*...
                        ( y(2).*w...
                         + thermal1d.yrToSec(obj.u)*obj.dTdx);
        end

        function dydx = thermalNonlinearSystem(obj, z, y, A_prev_func)
            % represent 1D heat transfer as a system of 1st order ODE
            % z is the independent variable (depth)
            % y(1) is T (temperature); y(2) is temperature gradient (dT/dz)
            switch obj.w_flag
                case "linear"
                    if obj.smb == 0; disp('SMB given = 0, hence w=0!');end
                    w = obj.calcLinearVerticalVel(z);
                case "DJ"
                    w = obj.calcDJVerticalVel(z);
                case "Nye"
                    w = obj.calcNyeVerticalVel(z);
                case "constant"
                    w = -thermal1d.yrToSec(obj.smb + obj.bm);
                otherwise
                    error("Vertical velocity request unknown!")
            end

            % Use the function handle to get A at the current z
            Q = obj.calcViscousDissipation(z, A_prev_func(z));
            
            dydx = zeros(2,1);
            dydx(1) = y(2);
            dydx(2) = y(2).*w./thermal1d.getAlpha()...
                      - Q./(thermal1d.rho_i*thermal1d.cp*thermal1d.getAlpha());
        end

        function [Pe, Peh, phi_hat] = non_dim_diffeq_coef(obj, w)
            %  Non-dimensionalize the variables
            % Arg:
            %    w: vertical velocity (z) 
            %    [not used yet for floating ice] Q: total heat source (shear heating + frictional heating at ice bed)
    
            Pe  = thermal1d.rho_i*thermal1d.cp*w*obj.thickness/thermal1d.k;
            Peh = thermal1d.rho_i*thermal1d.cp*obj.yrToSec(obj.u)...
                    *obj.thickness*(obj.thickness/thermal1d.L)/thermal1d.k;
            phi_hat = thermal1d.L*obj.dTdx/(obj.Ts-thermal1d.Tm);
            % Br = self.thickness**2 * Q/((self.Tm-self.Ts)*self.k)
            
        end
    
        function [theta_s, theta_b] = non_dim_bc(obj)
        % non-dimensionalize the boundary conditions
        % return:
        %    theta_s: non-dimensionalized ice surface temperature
        %    theta_b: non-dimensionalized ice base temperature (only if temperate)
        %    [not used yet for floating ice] dthetadz_b: non-dimensionalized basal temperature gradient (-GHF/k)
        
        % boundary conditions
        % F = self.calc_frictional_heat()
        
            theta_s = (obj.Ts-thermal1d.Tm)/(obj.Ts-thermal1d.Tm);
            theta_b = (thermal1d.Tm-thermal1d.Tm)/(obj.Ts-thermal1d.Tm);
            %dthetadz_b = -(self.G+F)*self.thickness/(self.k*(self.Ts-self.Tm))
        end

        function [theta, z_hat] = non_dim_diffeq_var(obj, T, z)
            % non-dimensionalize the variables in the diff eq.
            theta = (T-thermal1d.Tm)/(obj.Ts-thermal1d.Tm);
            z_hat = z/obj.thickness;
        end
        
        function [T, z] = dimensionalize(obj, theta, z_hat)
            % Recover the dimensions of temperature and axis
            T = theta*(obj.Ts-thermal1d.Tm) + thermal1d.Tm;
            z = z_hat*obj.thickness;
        end


        function [T,z] = solveMatrixFI(obj)
            % solve the floating ice shelf heat transfer 
            % using 2nd order finite difference by solving the linear sys.

            % here i am solving the non-dimensionalized version
            %     Pe(z)d(\theta)/dz - d^2 \theta /dz^2 + Pe_h = 0
            %     \theta(z=0) = 0 (base: at melting point)
            %     \theta(z=1) = 1 (surface)
            
            N = numel(obj.z);

            % get vertical velocity (it is the full profile/vector)
            % it is arranged from the surface to bottom
            switch obj.w_flag
                case "linear"
                    if obj.smb == 0; disp('SMB given = 0, hence w=0!');end
                    w = obj.calcLinearVerticalVel(obj.z);
                case "DJ"
                    w = obj.calcDJVerticalVel(obj.z);
                case "Nye"
                    w = obj.calcNyeVerticalVel(obj.z);
                case "constant"
                    w = -thermal1d.yrToSec(obj.smb + obj.bm)*ones(1,N);
                otherwise
                    error("Vertical velocity request unknown!")
            end
            
            % nondimensionalize the grid spacing
            h = obj.dz/obj.thickness;
            % allocate array space
            A = zeros(N,N);
            b = zeros(N,1);

            for ii = 2:N-1
                [Pe, Peh, phi_hat] = obj.non_dim_diffeq_coef(w(ii));
                if w(ii) > 0 % upwind for downward velocity
                   aa = -1/h^2;
                   bb = Pe/h + 2/h^2;
                   cc = -(Pe/h+1/h^2);

                else % w(ii) <= 0, upwind for upward velocity
                   aa = -1/h^2 + Pe/h;
                   bb = 2/h^2 - Pe/h;
                   cc = -1/h^2;
                end
                % assign
                A(ii,ii-1) = aa;
                A(ii,ii)   = bb;
                A(ii,ii+1) = cc;
                % RHS
                b(ii,1) = -Peh*phi_hat;
            end

            % apply boundary condition
            A(1,1) = 1;
            A(N,N) = 1;
            [theta_s, theta_b] = obj.non_dim_bc();
            b(1,1) = theta_s; 
            b(N,1) = theta_b;

            % solve 
            theta = linsolve(A,b);

            % bring back the dimension
            [T, z] = obj.dimensionalize(theta, obj.z/obj.thickness);

        end

        function init_g = guessInit(obj, z)
            % guess initial condition for both T and dT/dz
            % z is the position vector (not used in this simple guess)
            T_guess = 255; 
            gradT_guess = 30/1500; % 30 degree over 1.5 km 
            init_g = [T_guess; gradT_guess];
        end

        function obj = solveFI(obj)
            % Solve the 1D steady state linear advection-diffusion (BVP problem)
            % with bvp4c
            % FI stands for Floating Ice: base temperature at melting point
%             disp(' -----------------------')
%             disp('  Start boundary value problem solve for floating ice...')

            switch obj.solver
                case "bvp4c" % 4th order collocation method
                    % ice base is at the melting point
                    disp('    Solver: bvp4c')
                    solinit = bvpinit(obj.z, @(z) obj.guessInit(z));
                    obj.sol = bvp4c(@obj.thermalLinearSystemHAC,...
                                    @obj.applyBoundCondMeltBase,...
                                    solinit);
        
                    % disp('  Solution completed!')
                    % disp(' -----------------------')
                case "FD" % finite difference method
                    % disp('    Solver: 2nd order finite difference')
                    [obj.sol.y(1,:), obj.sol.x] = obj.solveMatrixFI();
                otherwise
                    error('Unknow solver')
            end
                    
        end
        
        function obj = solveLinearGI(obj)
            % Solve the 1D steady state linear advection-diffusion (BVP problem)
            % with bvp4c
            % GI stands for Grounded ice
            disp(' -----------------------')
            disp('  Start boundary value problem solve for grounded ice...')
            solinit = bvpinit(obj.z, @(z) obj.guessInit(z));
            obj.sol = bvp4c(@obj.thermalLinearSystemHAC,...
                            @obj.applyBoundCondSubmeltBase,...
                            solinit);

            if obj.sol.y(1,end)>thermal1d.Tm
                % the ice base should be at melting point
                % switch boundary condition
                disp('      Base exceeds melting point!')
                disp('      Entering new boundary value problem...')
                obj.sol = bvp4c(@obj.thermalLinearSystemHAC,...
                                @obj.applyBoundCondMeltBase,...
                                solinit);
            end

            disp('  Solution completed!')
            disp(' -----------------------')
        end

        function obj = solveNonlinearGI(obj)
            % Solve the 1D steady state nonlinear advection-diffusion (BVP problem)
            % with bvp4c; nonlinearity is addressed with Picard iteration
            % GI stands for Grounded Ice
            disp(' -----------------------')
            disp('  Start boundary value problem solve for grounded ice...')
            
            % Picard iteration parameters
            maxIter = 40;  % Maximum number of iterations
            tolerance = 1e-4;  % Convergence tolerance
            relaxation = 0.5;  % Relaxation factor
            
            % Initial guess (solve with constant A)
            solinit = bvpinit(obj.z, @(z) obj.guessInit(z));
            A_prev = thermal1d.computeGlenFlowRateParameter(260 * ones(size(obj.z)));

            % first make a linear solve
            obj.sol = bvp4c(@obj.thermalLinearSystem,...
                            @obj.applyBoundCondSubmeltBase,...
                            solinit);

            for i = 1:maxIter

                % create an interpolation function handle
                A_prev_func = @(z) interp1(obj.z, A_prev, z, 'linear','extrap');

                % Solve BVP with current A
                obj.sol = bvp4c(@(z,y) obj.thermalNonlinearSystem(z,y,A_prev_func), ...
                                @obj.applyBoundCondSubmeltBase, solinit);

                if obj.sol.y(1,end) > thermal1d.Tm

                    disp('      Base exceeds melting point!')
                    disp('      Entering new boundary value problem...')
                    obj.sol = bvp4c(@(z,y) obj.thermalNonlinearSystem(z,y,A_prev_func), ...
                                    @obj.applyBoundCondMeltBase, solinit);

                end

                % Extract temperature profile
                T_new = deval(obj.sol, obj.z);
                T_new(1,T_new(1,:)>thermal1d.Tm) = thermal1d.Tm; % no higher than pressure melting point

                % Calculate new A with relaxation
                A_new = (1-relaxation) * A_prev + relaxation * thermal1d.computeGlenFlowRateParameter(T_new(1,:));

                % Check for convergence
                rel_diff = norm(A_new - A_prev) / norm(A_prev);

                if mod(i, 5) == 0
                    disp(['  Picard iteration: ', num2str(i)])
                    disp(['      Relative difference in A: ', num2str(rel_diff)])
                end

                if rel_diff < tolerance
                    disp(['  Converged after ', num2str(i), ' iterations'])
                    break;
                end

                A_prev = A_new;

                solinit = obj.sol;  % Use current solution as initial guess for next iteration
                if i == maxIter
                    warning('Maximum iterations reached without convergence')
                end
            end


            disp('  Solution completed!')
            disp(' -----------------------')
        end

                
        function w = calcLinearVerticalVel(obj, z)
            % Get the vertical velocity, assuming linear in depth and 0
            % basal melting
            w = -obj.smb.*z./obj.thickness;
            w = thermal1d.yrToSec(w);
        end

        function w = calcDJVerticalVel(obj, z)
            % get the vertical velocity, assuming Dansgaard-Johnsen
            % formulation, which is applicable to ice divide
            h = 0.5*obj.thickness;
            a_norm = 2*obj.smb/(2*obj.thickness - h); % a_norm is called k in the original formulation
            
            w = zeros(size(z));
            lower_idx = z>=0 & z<=h;
            upper_idx = z<=obj.thickness & z > h;
            w(lower_idx) = -a_norm*z(lower_idx).^2/(2*h) - obj.bm*(1-z(lower_idx)/obj.thickness);
            w(upper_idx) = -a_norm/2*(2*z(upper_idx) - h) - obj.bm*(1-z(upper_idx)/obj.thickness);
            
            % unit conversion
            w = thermal1d.yrToSec(w);
        end

        function w = calcNyeVerticalVel(obj, z)
            % get the vertical velocity, assuming Nye model, which is a
            % linear model but accounting for basal melt rate
            w = -(z*obj.smb/obj.thickness + (1 - z/obj.thickness)*obj.bm);

            % unit conversion
            w = thermal1d.yrToSec(w);
        end
            
        
        function Q = calcViscousDissipationFromT(obj, z, T)
            % calculate the heat source from viscous dissipation
            % prefactor A is computed from temperature
            A = thermal1d.computeGlenFlowRateParameter(T);
            dudz = 2*A.*(thermal1d.rho_i*thermal1d.g*abs(obj.dsdx).*(obj.thickness - z)).^thermal1d.n;
            Q = 2*((dudz).^(thermal1d.n+1)./A).^(1./thermal1d.n);

        end

        function Q = calcViscousDissipation(obj, z, A)
            % Get the heat source at depths from viscous dissipation
            % ...provided A (used for Picard iteration solving the
            % nonlinear version)
            dudz = 2*A.*(thermal1d.rho_i*thermal1d.g*abs(obj.dsdx).*(obj.thickness - z)).^thermal1d.n;
            Q = 2*((dudz).^(thermal1d.n+1)./A).^(1./thermal1d.n);
%             disp(['Q term1:' num2str((dudz).^(thermal1d.n+1)./A), 'A:', num2str(A)])
        end

        
        function bc_res = applyBoundCondMeltBase(obj, ya, yb)
            % Apply a boundary condition (works with MATLAB bvp4c)
            % Assuming the base is at melting point (i.e. Dirichlet
            % condtions at both the surface and base)
            %
            % left boundary: ice surface, z = H
            % right boundary: ice base, z = 0
            bc_res = [ya(1) - obj.Ts
                      yb(1) - thermal1d.Tm];
        end

        function bc_res = applyBoundCondSubmeltBase(obj, ya, yb)
            % Apply a boundary condition (works with MATLAB bvp4c)
            % Assuming the base is below the melting point
            % hence at the ice base we apply geothermal heat flux (Neumann
            % condition)

            % left boundary: ice surface, z = H
            % right boundary: ice base, z = 0
            bc_res = [ya(1) - obj.Ts
                      yb(2) + obj.G/thermal1d.k];
        end


        function printParameters(obj)
            % Print the parameters of the thermal1d model
            fprintf('Thermal1D Model Parameters:\n')
            fprintf('----------------------------\n')
            fprintf('Ice sheet thickness: %.2f m\n', obj.thickness)
            fprintf('Surface slope: %.4f\n', obj.dsdx)
            fprintf('Vertical resolution (dz): %.2f m\n', obj.dz)
            fprintf('Surface mass balance: %.2f m/year\n', obj.smb)
            fprintf('Surface temperature: %.2f K\n', obj.Ts)
            fprintf('Basal mass balance: %.2f m/year\n', obj.bm)
            fprintf('Geothermal heat flux: %.4f W/m^2\n', obj.G)
            fprintf('\nConstant Parameters:\n')
            fprintf('Thermal conductivity (k): %.2f W/(m·K)\n', thermal1d.k)
            fprintf('Specific heat capacity (cp): %.2f J/(kg·K)\n', thermal1d.cp)
            fprintf('Ice density (rho_i): %.2f kg/m^3\n', thermal1d.rho_i)
            fprintf('Glen''s flow law exponent (n): %d\n', thermal1d.n)
            fprintf('Gravitational acceleration (g): %.2f m/s^2\n', thermal1d.g)
            fprintf('Melting point (Tm): %.2f K\n', thermal1d.Tm)
            fprintf('Thermal diffusivity (alpha): %.2e m^2/s\n', thermal1d.getAlpha())
        end

    
        function plot(obj)
            % Plot the temperature profile

            plot(obj.sol.y(1,:), obj.sol.x,'LineWidth',2)
            xlabel('Temperature (K)')
            ylabel('Elevation (m)')
        end

    end

    methods (Static)
        function A = computeGlenFlowRateParameter(T)
            % Compute the flow law parameter A
            % unit: Pa^-3 s^-1
            [~,Qc,~] = thermal1d.defineActivationEnergies(T);
            R = 8.314; % J/mol K
            A0 = 2.4e-24./exp(-(115000./R).*((1/273)-(1/263)));

            A = zeros(size(T));
            if length(A0)>1
                for i=1:length(T)
                    A(i) = A0(i).*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
                end
            else
                for i=1:length(T)
                    A(i) = A0.*exp(-(Qc(i)./R).*((1/T(i))-(1/263)));
                end
            end
        end

        function alpha = getAlpha()
            % Calculate and return alpha
            alpha = thermal1d.k / (thermal1d.rho_i * thermal1d.cp);
        end

        function var = yrToSec(var)
            % convert yr to sec
            var = var/(3.154e7);
        end

        function [Qg,Qc,Qm] = defineActivationEnergies(T)
        
            Temp = T-273;
            
            % Activation Energy for Creep
            tp = 0;
            tm = -20;
            tc = -10;
            Qcp = 100; %kJ/mol
            Qcm = 60; %kJ/mol
            c1 = (Qcp-Qcm)./(atan(tp-tc)-atan(tm-tc));
            c2 = Qcp-c1*tanh(tp-tc);
            
            % Activation Energy for Grain Growth
            % define Arrhenius relation for grain growth
            tp = 0;
            tm = -20;
            tc = -10;
            Qgp = 100; %kJ/mol
            Qgm = 40; %kJ/mol
            g1 = (Qgp-Qgm)./(atan(tp-tc)-atan(tm-tc));
            g2 = Qgp-g1*tanh(tp-tc);
            
            % Activation Energy for Grain Boundary Mobility
            % define Arrhenius relation for grain boundary mobility
            tp = 0;
            tm = -20;
            tc = -10;
            Qmp = 40; %kJ/mol
            Qmm = 100; %kJ/mol
            m1 = (Qmp-Qmm)./(atan(tp-tc)-atan(tm-tc));
            m2 = Qmp-m1*atan(tp-tc);
            
            Qg = zeros(size(T));
            Qg = g1*atan(Temp-tc)+g2;
            Qg = Qg.*1e3;
            
            Qc = zeros(size(T));
            Qc = c1*atan(Temp-tc)+c2;
            Qc = Qc.*1e3;
            
            % define Arrhenius relation for grain boundary mobility
            Qm = zeros(size(T));
            Qm = m1*atan(Temp-tc)+m2;
            Qm = Qm.*1e3;
        
        end

    end
end
