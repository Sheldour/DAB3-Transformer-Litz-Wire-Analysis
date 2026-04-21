% ---------------------------------------------------------
% This script implements the analytical model for calculating the optimum Litz-wire strand count, as proposed in the following paper:
% [1] C. R. Sullivan, “Optimal choice for number of strands in a litz-wire transformer winding,” IEEE Transactions on Power Electronics, vol. 14, no. 2, pp. 283–291, Mar. 1999, doi: 10.1109/63.750181.

% ---------------------------------------------------------
% Compute optimum Litz-wire strand count using Sullivan's closed-form model, parameters from the paper's example in section V.

clear; clc;

%% --------------------- User Inputs ---------------------
% Core / transformer (prototype references)
f_s = 375e3;               % Switching frequency [Hz]
omega = 2 * pi * f_s;     % Angular frequency [rad/s]
N = 14;                   % Turns (primary)

% EEW58 window data
b_c = 6.3e-3;              % Core window breadth B-C [m] (EEW58)
b_b = 4.93e-3;        % Bobbin breadth for the considered winding [m]
h = 1.09e-3;      % Height allocated to this winding [m]

% Packing factors
F_p = 0.85;               % Basic packing factor (set by winding layout)
F_lp = 0.505;              % Litz-process factor (serving/bundle/twist effects)

% Conductor and model constants
mu_0 = 4 * pi * 1e-7;     % Vacuum permeability [H/m]
rho_c = 1.68e-8;         % Copper resistivity [ohm*m]
k = 1.0;                  % Field-distribution factor; usually 1
d_r = 0.079e-3;           % Reference diameter [m] (AWG40, 0.079 mm)

% Skin depth at switching frequency
delta_skin = sqrt(rho_c / (pi * mu_0 * f_s));   % [m]

% Insulation type constants (from paper text)
insulation_type = "single";  % "single" or "heavy"
switch lower(insulation_type)
    case "single"
        beta = 0.97;
        alpha = 1.12;
    case "heavy"
        beta = 0.94;
        alpha = 1.24;
    otherwise
        error('Unknown insulation_type. Use "single" or "heavy".');
end

%% --------------------- Core Equations ---------------------
% NOTE:
% In some OCR copies, denominator appears as p^2. Based on symbol list and
% original model consistency, rho_c^2 is used here.

gamma_val = (pi^2 * N^2 * omega^2 * mu_0^2 * d_r^(6 - 6 / beta) * alpha^(-6 / beta) ...
    * (F_p * F_lp * b_b * h / N)^(3 / beta) * k) ...
    / (768 * rho_c^2 * b_c^2);

n_opt = (((2 / beta - 1) * gamma_val) / (1 / beta - 1))^(1 / (3 / beta - 2));

if ~isreal(n_opt) || n_opt <= 0
    error('Computed n_opt is non-physical. Check inputs and units.');
end

% Diameter formula from the same section snapshot:
% d_t = sqrt(F_p * F_lp * b_b * h / (n * N))
d_t_opt = sqrt(F_p * F_lp * b_b * h / (n_opt * N));

if d_t_opt * 1e3 < 0.00786 || d_t_opt * 1e3 > 0.255
    warning(['d_t at n_opt is outside the 30-60 AWG validity range ', ...
        '(0.00786-0.255 mm). Check area allocation and packing factors.']);
end

% Practical strand-count suggestion
preferred_counts = [20 25 30 40 50 60 70 80 90 100 110 120 125 130 140 150 ...
    160 180 200 220 240 250 280 300 320 350 400 450 500 560 600 700 800 900 1000];
[~, idx] = min(abs(preferred_counts - n_opt));
n_preferred = preferred_counts(idx);

%% --------------------- Sensitivity to F_lp ---------------------
F_lp_vec = 0.5:0.05:0.8;
n_vec = zeros(size(F_lp_vec));
for i = 1:numel(F_lp_vec)
    Fi = F_lp_vec(i);
    gamma_i = (pi^2 * N^2 * omega^2 * mu_0^2 * d_r^(6 - 6 / beta) * alpha^(-6 / beta) ...
        * (F_p * Fi * b_b * h / N)^(3 / beta) * k) ...
        / (768 * rho_c^2 * b_c^2);
    n_vec(i) = (((2 / beta - 1) * gamma_i) / (1 / beta - 1))^(1 / (3 / beta - 2));
end

%% --------------------- Print Results ---------------------
fprintf('===== Optimum Litz Strand Count (Sullivan model) =====\n');
fprintf('insulation_type = %s\n', insulation_type);
fprintf('alpha = %.4f, beta = %.4f\n', alpha, beta);
fprintf('f_s = %.1f kHz, N = %d\n', f_s / 1e3, N);
fprintf('b_c (EEW58, B-C) = %.2f mm\n', b_c * 1e3);
fprintf('b_b (bobbin breadth) = %.2f mm\n', b_b * 1e3);
fprintf('h (allocated to this winding) = %.2f mm\n', h * 1e3);
fprintf('F_p = %.3f, F_lp = %.3f\n', F_p, F_lp);
fprintf('skin depth delta = %.4f mm\n', delta_skin * 1e3);
fprintf('gamma = %.6e\n', gamma_val);
fprintf('n_opt (real) = %.3f\n', n_opt);
fprintf('n_opt (nearest integer) = %d\n', round(n_opt));
fprintf('n_opt (preferred count) = %d\n', n_preferred);
fprintf('d_t at n_opt = %.4f mm\n', d_t_opt * 1e3);

fprintf('\nF_lp sweep (for quick check):\n');
for i = 1:numel(F_lp_vec)
    fprintf('  F_lp = %.2f -> n_opt = %.2f\n', F_lp_vec(i), n_vec(i));
end
