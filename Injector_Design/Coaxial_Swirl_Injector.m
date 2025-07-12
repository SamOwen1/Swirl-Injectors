% -------------------------------------------------------------------------
% Coaxial Swirl Injector Design Script %
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

% Choose Coefficients of Nozzle Opening %
% Open Stages (Eta<1) - Rv = Rn
% Closed Stages (Eta>1) - Rv > Rn
eta1 = 3; % Inner Element Coefficient of Nozzle Opening
eta2 = 1.2; % Outer Element Coefficient of Nozzle Opening

% Fluid Parameters %
mdot1 = 0.04; % Mass Flow Rate of Inner Element (kg/s)
mdot2 = 0.2; % Mass Flow Rate of Outer Element (kg/s)
deltap1 = 8e5; % Pressure Drop of Inner Element (Pa)
deltap2 = 5e5; % Pressure Drop of Outer Element (Pa)
rho1 = 1000; % Fluid Density of Inner Element (kg/m^3)
rho2 = 1000; % Fluid Density of Outer Element (kg/m^3)
n1 = 3; % Number of Tangential Inlets of Inner Element
n2 = 5; % Number of Tangential Inlets of Outer Element 
alpha1 = 110; % Full Spray Cone Angle of Inner Element (degrees)
xi1 = 0.5; % Hydraulic Loss Coefficient of Inner Element - Set to 0 for No Loss
xi2 = 0.5;

op_cl1 = 0; % 0 for Closed Inner Element 1 for Open Inner Element

if op_cl1 == 0

    % Calculate Inner Element Filling Efficiency %
    phi_eq1 = @(phi1) ((sqrt(2) * (1 - phi1)) / sqrt(2 - phi1)) - sin(alpha1 * pi / 360);
    phi1 = fzero(phi_eq1, 0.8);

    % Determine Inner Element Flow Coefficient %
    A1 = (sqrt(2) * (1 - phi1)) / (phi1 * sqrt(phi1));
    mu1 = sqrt(((phi1 ^ 3) * (eta1 ^ 2)) / ((2 - phi1) * (eta1 ^ 2) + (xi1 * A1 ^ 2 * phi1 ^ 3))); 

    % Inner Element Geometry %
    rn1 = sqrt(mdot1 / (pi * mu1 * sqrt(2 * rho1 * deltap1)));
    rin1 = rn1 * eta1;
    rt1 = sqrt((rn1 * rin1) / (n1 * A1));
    rv1 = rin1 + rt1;

elseif op_cl1 == 1

    A1 = eta1 / (n1 * (1 - eta1) ^ 2); % A Which Produces Case rv1=rn1
    phi1 = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A1, 0.8); 
    mu1 = sqrt(((phi1 ^ 3) * (eta1 ^ 2)) / ((2 - phi1) * (eta1 ^ 2) + (xi1 * A1 ^ 2 * phi1 ^ 3)));
    alpha1 = (360 / pi) * asin(((sqrt(2) * (1 - phi1)) / sqrt((2 - phi1))));

    % Inner Element Geometry %
    rn1 = sqrt(mdot1 / (pi * mu1 * sqrt(2 * rho1 * deltap1)));
    rin1 = rn1 * eta1;
    rt1 = sqrt((rn1 * rin1) / (n1 * A1));
    rv1 = rin1 + rt1;

end

% Inner Element Outer Wall Radius %
delta = 0.5e-3; % Wall Thickness 0.2mm<delta<0.8mm (m)
r1o = rn1 + delta;

% Outer Element Design %
op_cl2 = 0; % 0 for Closed Outer Element 1 for Open Outer Element

if op_cl2 == 0

    d = 10; % Usually 10<d<15
    alpha2 = 0.5 * ((2 * alpha1) - d); % Full Spray Cone Angle of Outer Element (degrees)

    % Calculate Outer Element Filling Efficiency %
    phi_eq2 = @(phi2) ((sqrt(2) * (1 - phi2)) / sqrt(2 - phi2)) - sin(alpha2 * pi / 360);
    phi2 = fzero(phi_eq2, 0.8);

    % Determine Outer Element Flow Coefficient %
    A2 = (sqrt(2) * (1 - phi2)) / (phi2 * sqrt(phi2));
    mu2 = sqrt(((phi2 ^ 3) * (eta2 ^ 2)) / ((2 - phi2) * (eta2 ^ 2) + (xi2 * A2 ^ 2 * phi2 ^ 3))); 

    % Outer Element Geometry %
    rn2 = sqrt(mdot2 / (pi * mu2 * sqrt(2 * rho2 * deltap2)));
    rin2 = rn2 * eta2;
    rt2 = sqrt((rn2 * rin2) / (n2 * A2));
    rv2 = rin2 + rt2;

elseif op_cl2 == 1

    A2 = eta2 / (n2 * (1 - eta2) ^ 2); % A Which Produces Case rv2=rn2
    phi2 = fzero(@(phi) ((sqrt(2) * (1 - phi)) / (phi * sqrt(phi))) - A2, 0.8); 
    mu2 = sqrt(((phi2 ^ 3) * (eta2 ^ 2)) / ((2 - phi2) * (eta2 ^ 2) + (xi2 * A2 ^ 2 * phi2 ^ 3)));
    alpha2 = (360 / pi) * asin(((sqrt(2) * (1 - phi2)) / sqrt((2 - phi2))));

    % Outer Element Geometry %
    rn2 = sqrt(mdot2 / (pi * mu2 * sqrt(2 * rho2 * deltap2)));
    rin2 = rn2 * eta2;
    rt2 = sqrt((rn2 * rin2) / (n2 * A2));
    rv2 = rin2 + rt2;

end

% Find Fluid Radii of Outer Element %
rn_l2 = rn2 * sqrt(1 - phi2); % Nozzle 
phiv0 = ((3 * phi2) - (2 * phi2 ^ 2)) / (2 - phi2);
rv0_l2 = rn2 * sqrt(1 - phiv0); % Head End
rv0_l2_norm = sqrt(1 - phiv0);
rv2_norm = rv2 / rn2;
rv_l2 = rn2 * fzero(@(rv) ((rv0_l2_norm^2 * (rv2_norm^2 - rv^2)) / (rv2_norm^2 - rv^2 - mu2^2)) - (rv^2), 0.5); % Vortex Chamber

% Check Gas Vortex Radius - Assumes Vortex Chamber of Outer Element is Adjacent to Inner Element Nozzle
% Change Geometry Checks for Different Layouts
if rn_l2 - r1o < 0

    gvn = rn_l2 - r1o;
    fprintf('Inner Element Not Accommodated: %.4f m', gvn);
    fprintf('DESIGN INVALID - Change Parameters');

elseif rv0_l2 - r1o < 0

    gvv0 = rv0_l2 - r1o;
    fprintf('Inner Element Not Accommodated: %.4f m', gvv0);
    fprintf('DESIGN INVALID - Change Parameters');

elseif rv_l2 - r1o < 0

    gvv = rv_l2 - r1o;
    fprintf('Inner Element Not Accommodated: %.4f m', gvv);
    fprintf('DESIGN INVALID - Change Parameters');

else

    fprintf('DESIGN IS VALID')
    % Aim for 0.3mm
    gvn = rn_l2 - r1o;
    fprintf('Nozzle Gas Vortex: %.4f mm', gvn * 1e3);
    gvv0 = rv0_l2 - r1o;
    fprintf('Head End Gas Vortex: %.4f mm', gvv0 * 1e3);
    gvv = rv_l2 - r1o;
    fprintf('Vortex Chamber Gas Vortex: %.4f mm', gvv * 1e3);

end

Km = mdot2 / mdot1; % Propellent Ratio
tao = 0.2e-3; % Mixing Time 0.2ms for Non-Hypergolic - 0.1ms for Hypergolic
lmix = sqrt(2)*tao*(((Km*mu2/((Km+1)*phi2)) * sqrt(deltap2/rho2)) +  ((mu1/((Km+1)*phi1)) * sqrt(deltap1/rho1)));
del = ((rn2 - rn1) / tan(alpha1 * pi /360));

rec = 0; % 0 For External Mixing 1 For Impinged 2 for Internal Mixing

if rec == 0

    recess = del;
    fprintf('0 mm < Recess < %.4f mm', recess * 1e3);

elseif rec == 1

    recess = del;
    fprintf('Recess = %.4f mm', recess * 1e3);

elseif rec == 2
   
    recess = del + lmix;
    fprintf('%.4f mm < Recess < %.4f mm', del * 1e3, recess * 1e3);

end

% Final Checks

if rv1 < rn1

    fprintf('DESIGN INVALID - Change Parameters');

elseif rv2 < rn2

    fprintf('DESIGN INVALID - Change Parameters');

elseif rec == 0 && alpha2 >= alpha1

    fprintf('DESIGN INVALID - Change Parameters');

elseif op_cl2 == 1 && rv2 > rn2

    fprintf('DESIGN INVALID - Change Parameters');

elseif op_cl1 == 1 && rv1 > rn1

    fprintf('DESIGN INVALID - Change Parameters');

end

% Print Geometry
fprintf('INNER ELEMENT:')
fprintf('Nozzle Radius %.4f mm', rn1 * 1e3);
fprintf('Vortex Chamber Radius %.4f mm', rv1 * 1e3);
fprintf('Tangential Inlet Radius %.4f mm', rt1 * 1e3);
fprintf('Outer Wall Radius (Nozzle) %.4f mm', r1o * 1e3);

fprintf('OUTER ELEMENT:')
fprintf('Nozzle Radius %.4f mm', rn2 * 1e3);
fprintf('Vortex Chamber Radius %.4f mm', rv2 * 1e3);
fprintf('Tangential Inlet Radius %.4f mm', rt2 * 1e3);
