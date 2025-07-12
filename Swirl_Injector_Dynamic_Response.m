% -------------------------------------------------------------------------
% Bazarov Swirl Injector Response %
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

% Define Given Parameters %
Rn = 0.002755978; % Nozzle Radius (m)
Rv = 3 * Rn; % Vortex Chamber Radius (m)
Rt = 0.460069 * Rn; % Tangential Inlet Radius (m)
lv = 25e-3; % Uniform Vortex Chamber Length (m)
ln = 8e-3; % Nozzle Length (m)
lc = 0e-3; % Non-Uniform Vortex Chamber Length (m)
lt = 2e-3; % Tangential Inlet Length (m)
n = 4; % Number of Tangential Inlets
rho = 1000; % Fluid density (kg/m^3)
deltap = 0.6e6; % Pressure Drop (Pa)

% Derived Variables %
Rin = Rv - Rt; 
eta = Rin / Rn; % Coefficient of Nozzle Opening
A = (Rn * Rin) / (n * Rt ^ 2); % Geometric Constant
phi_guess = 0.8;
phi = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A, phi_guess); % Filling Efficiency
mu = sqrt((phi ^ 3) / (2 - phi)); % Mass Flow Coefficient
usum = sqrt((2 * deltap) / rho); % Total Velocity
phiv0 = ((3 * phi) - (2 * phi ^ 2)) / (2 - phi); % Head End Filling Efficiency
rn_norm = sqrt(1 - phi); % Normalized Nozzle Film Radius
rv0_norm = sqrt(1 - phiv0); % Normalized Head End Film Radius
Rn_norm = Rn / Rn; % Normalized Nozzle Radius
Rv_norm = Rv / Rn; % Normalized Vortex Chamber Radius

% Solve Normalized Vortex Chamber Film Radius %
rv_norm = fzero(@(rv) ((rv0_norm^2 * (Rv_norm^2 - rv^2)) / (Rv_norm^2 - rv^2 - mu^2)) - (rv^2), phi_guess);

% Solve Axial Velocities %
uvz = usum * sqrt(1 - (rv0_norm ^ 2) / (rv_norm ^ 2));
unz = usum * sqrt(1 - (rv0_norm ^ 2) / (rn_norm ^ 2));
uin = usum * rv0_norm / eta;

% Angular Momentum Constant %
C = usum * rv0_norm * Rn; 

% Define Normalized Variables %
Rt_norm = Rt / Rn;
ln_norm = ln / Rn;
lv_norm = lv / Rn;
lt_norm = lt / Rn;
lc_norm = lc / Rn;
C_norm = C / (Rn * uin);
uvz_norm = uvz / uin;
unz_norm = unz / uin;
uin_norm = uin / uin;

% Wave Speeds %
cv_norm = uvz_norm + sqrt((C_norm ^ 2) * ((Rv_norm ^ 2) - (rv_norm ^ 2)) / (2 * rv_norm ^ 4));
cn_norm = unz_norm + sqrt((C_norm ^ 2) * ((Rn_norm ^ 2) - (rn_norm ^ 2)) / (2 * rn_norm ^ 4));

% Define Frequency Range (Hz) %
freq_range = linspace(0, 2200, 2200);
omega_range = freq_range .* (2 * pi); % Convert to Angular Frequency
omega_range_norm = (Rn .* omega_range) ./ uin;
response_magnitude = zeros(size(freq_range)); % Preallocate Respones Magnitude
nu = 0.1; % Artificial Viscosity Coefficient

% Loop Over Frequencies %
for idx = 1:length(omega_range_norm)
    omega = omega_range_norm(idx);

    % Compute Phase Shifts %
    phiv = omega * ((lv_norm + lc_norm) / cv_norm);
    phin = omega * (ln_norm / cn_norm);

    % Compute Reflection Coefficient %
    Pi_refl = 1 - 2 * (sqrt(phi) / sqrt((Rv_norm ^ 2) - (rv0_norm ^ 2)));

    % Compute Surface Wave Transfer Functions Using Infinite Sum Formulas %
    V2_constant = 1 / (A * sqrt(2 * (Rv_norm ^ 2 - rv0_norm ^ 2)));
    sum_terms_V2 = 1 / (1 - Pi_refl * exp(-2 * phiv * (1i + nu)));
    sum_terms_VN = exp(-phiv * (1i + nu)) / (1 - Pi_refl * exp(-2 * phiv * (1i + nu)));

    % Compute Pi_V2 and Pi_VN %
    Pi_V2 = V2_constant * sum_terms_V2;
    Pi_VN = sum_terms_VN;

    % Define Strouhal Numbers %
    Sh_V = omega * Rv_norm / uvz_norm; 
    Sh_T = omega * lt_norm / uin_norm; 

    % Define Function for Vorticity Wave Transfer Function %
    f_x = @(x) x .* Sh_V .* tan(pi * x / 2);
    K = 1 - (rv0_norm / Rv_norm);

    % Compute Real and Imaginary Parts of Pi_V3 Using Numerical Integration %
    Re_Pi_V3 = 2 * integral(@(x) (cos(f_x(x)) ./ (1 - K * x).^3) .* exp(-nu * f_x(x)), 0, 1, 'ArrayValued', true);
    Im_Pi_V3 = -2 * integral(@(x) (sin(f_x(x)) ./ (1 - K * x).^3) .* exp(-nu * f_x(x)), 0, 1, 'ArrayValued', true);
    Pi_V3 = Re_Pi_V3 + 1i * Im_Pi_V3;

    % Compute Nozzle Transfer Function %
    Pi_N = (1 - Pi_refl) * exp(-1i * phin);

    % Compute Tangential Inlet Responce Function %
    Pi_T = 0.5 * (1 - 1i * Sh_T) / (1 + Sh_T^2);

    % Compute Total Injector Responce Function %
    Pi_Sigma = (Rv_norm / rv0_norm)^2 * (Pi_T * Pi_VN * Pi_N) / (1 + 2 * Pi_T * (Pi_V2 + Pi_V3));
    response_magnitude(idx) = abs(Pi_Sigma);
end

% Plotting %
figure;
plot(freq_range, response_magnitude, 'b', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
xlim([0 2200])
ylabel('|Π_{\Sigma}|');
grid on;
hold off;