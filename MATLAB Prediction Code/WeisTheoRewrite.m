%% Rewrite of the Weisshaar Metrix method for the theodorsen function
% Alexander Ketzle
clc, clear
i = sqrt(-1);

% page 219 parameters
% mu = 76;
% a_h = -0.15;
% x_alpha = 0.25;
% r_alpha = sqrt(0.388);
% b = 5 / 12; % ft
% freq_alpha = 64.1; % rad/s
% freq_h = 55.9; % rad/s
% g_alpha = 0;
% g_h = 0;

% doppler's fin
% a_h =  0.0;
% b = 2.25/12;
% x_alpha = 0.0;
% r_alpha = 0.57785;
% freq_alpha = 2866.77844;
% freq_h = 1853.69007;
% mu = 110.28942;
% g_alpha = 0.005;
% g_h = 0.005;

% Cippola N5800 Parameters - This is one of the examples shipped with
%FinSim
a_h = 0.0;
b = 3.5625 / 12; % ft
x_alpha = 0.316;
r_alpha = 0.57757;
freq_alpha = 2593.373; % rad/s
freq_h = 2458.08525; % rad/s
mu = 77.11441;
g_alpha = 0.00;
g_h = 0.00;

%% Calculation

kStepSize = 0.0001;
kRange = [kStepSize,14];
n = ((kRange(2) - kRange(1)) / kStepSize) + 1;
kSet = linspace(kRange(1),kRange(2),n);

k = kSet;

Ch_k = besselh(1,2,k) ./ (besselh(1,2,k) + (i * besselh(0,2,k))); % Theodorsen Function; Less lines to compute than using the other bessel functions
F = real(Ch_k);
G = imag(Ch_k);

L_h = 1 - (i .* 2 .* Ch_k ./ k);
L_alpha = 0.5 - (i .* (1 + (2 .* Ch_k)) ./ k) - (2 .* Ch_k ./ k.^2);
M_alpha = (3/8) - (i ./ k);

A11 = (freq_alpha.^2 ./ freq_h.^2) .* (1 + (L_h ./ mu));
A12 = (freq_alpha.^2 ./ freq_h.^2) .* (x_alpha + (mu.^-1 .* (L_alpha - ((0.5 + a_h) .* L_h))));
A21 = ((x_alpha ./ r_alpha.^2) + ((1 ./ (mu .* r_alpha.^2)) .* (0.5 - (0.5 + a_h) .* L_h)));
A22 = 1 + ((1 ./ (mu .* r_alpha.^2)) .* (M_alpha - ((L_alpha + 0.5) .* (0.5 + a_h)) + (L_h .* (0.5 + a_h).^2)));

fluttermatrix = [A11, A12; A21, A22];
omega = eig(fluttermatrix);