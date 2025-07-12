% -------------------------------------------------------------------------
% Single Element Swirl Injector Code
% Author: Sam Owen
% License: MIT License
% -------------------------------------------------------------------------

% Copyright (c) 2025 Sam Owen

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE

clc;
clear;

% Define Variables
mdot = 0.533646 / 15; % Injector Mass Flow Rate (kg/s)
rho = 786; % Fluid Density (kg/m^3)
deltap = 5e5; % Total Pressure Drop Across Injector (Pa)
n = 3; % Number of Tangential Inlets

% Define Solution Arrays
phis = linspace(0.1, 0.9, 100); % Filling Efficiency
eta_values = [0.7, 1.4, 2.1]; % Coefficient of Nozzle Opening
xi = 0.5; % Hydraulic Loss Coefficient for Sharp Edged Inlet - Set to 0 for No Loss

rns_all = zeros(length(phis), length(eta_values)); % Nozzle Radius Array
rvs_all = zeros(length(phis), length(eta_values)); % Vortex Chamber Radius Array
rts_all = zeros(length(phis), length(eta_values)); % Tangential Inlet Radius Array
alphas_epsilon_all = zeros(length(phis), length(eta_values)); 

% Perform Calculations
for eta_idx = 1:length(eta_values)
    eta = eta_values(eta_idx);
    
    for i = 1:length(phis)
        phi = phis(i);

        % Geometric Constant % 
        A = (sqrt(2) * (1 - phi)) / (phi * sqrt(phi));

        % Hydraulic Mass Flow Coefficient %
        mu = sqrt(((phi ^ 3) * (eta ^ 2)) / ((2 - phi) * (eta ^ 2) + (xi * A ^ 2 * phi ^ 3))); 
        
        % Nozzle Radius %
        rn = sqrt(mdot / (pi * mu * sqrt(2 * rho * deltap)));
        rin = rn * eta;

        % Tangential Inlet Radius %
        rt = sqrt((rn * rin) / (n * A));

        % Vortex Chamber Radius % 
        rv = rin + rt;
        
        % Hydraulic Total Velocity %
        usum = sqrt((2 * deltap) / rho) * sqrt(1 - ((xi * mu ^ 2 * A ^ 2) / (eta ^ 2))); 

        % No Loss Mass Flow Coefficient %
        mu0 = sqrt((phi ^ 3) / (2 - phi));

        % Exit Spray Half Angle %
        alphae = atan((mu0 * A) / sqrt(1 - ((mu0 ^ 2) * (A ^ 2))));
        
        % Store Results %
        rns_all(i, eta_idx) = rn;
        rvs_all(i, eta_idx) = rv;
        rts_all(i, eta_idx) = rt;
        alphas_epsilon_all(i, eta_idx) = alphae;

    end
end

% Plot Results %
figure;
hold on;
for eta_idx = 1:length(eta_values)
    plot(phis, rns_all(:, eta_idx) * 1e3, 'LineWidth', 2);
end
ylabel('r_N (mm)');
xlabel('φ');
legend(arrayfun(@(eta) sprintf('η = %.1f', eta), eta_values, 'UniformOutput', false));
grid on;
hold off;

figure;
hold on;
for eta_idx = 1:length(eta_values)
    plot(phis, rvs_all(:, eta_idx) * 1e3, 'LineWidth', 2);
end
ylabel('r_V (mm)');
xlabel('φ');
legend(arrayfun(@(eta) sprintf('η = %.1f', eta), eta_values, 'UniformOutput', false));
grid on;
hold off;

figure;
hold on;
for eta_idx = 1:length(eta_values)
    plot(phis, rts_all(:, eta_idx) * 1e3, 'LineWidth', 2);
end
ylabel('r_T (mm)');
xlabel('φ');
legend(arrayfun(@(eta) sprintf('η = %.1f', eta), eta_values, 'UniformOutput', false));
grid on;
hold off;

figure;
hold on;
for eta_idx = 1:length(eta_values)
    plot(phis, 2 * alphas_epsilon_all(:, eta_idx) * (180 / pi), 'LineWidth', 2);
end
ylabel('2\alpha (degrees)');
xlabel('φ');
legend(arrayfun(@(eta) sprintf('η = %.1f', eta), eta_values, 'UniformOutput', false));
grid on;
hold off;
