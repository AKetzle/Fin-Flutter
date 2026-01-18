%% Y.C. Fung Ch. 6.9 Flutter Calculation
% based on Theodorsen Method
% Alexander Ketzle
% starts around page 213
clc, clear, close all;
i = sqrt(-1);

% page 219 parameters
mu = 76;
a_h = -0.15;
x_alpha = 0.25;
r_alpha = sqrt(0.388);
b = 5 / 12; % ft
freq_alpha = 64.1; % rad/s
freq_h = 55.9; % rad/s
g_alpha = 0;
g_h = 0;

k = 0;
kstep = 0.01;
krange = [kstep,6];


Ch_k = besselh(1,2,k) / (besselh(1,2,k) + (i * besselh(0,2,k))); % Theodorsen Function; Less lines to compute than using the other bessel functions
F = real(Ch_k);
G = imag(Ch_k);

A_R = -(mu + 1) - (2 * G / k);
A_I = 2 * F / k;
B_R = -((mu * x_alpha) - a_h) + (2 * F / k^2) - ((0.5 - a_h) * 2 * G / k);
B_I = (1 / k) * (1 + (2 * G / k) + ((0.5 - a_h) * 2 * F));
D_R = -((mu * x_alpha) - a_h) + ((0.5 + a_h) * 2 * G / k);
D_I = -(0.5 + a_h) * 2 * F / k;
E_R = -((mu * r_alpha^2) + a_h^2 + 0.125) + ((0.25 - a_h^2) * 2 * G / k) - ((0.5 + a_h) * 2 * F / k^2);
E_I = (1 / k) * ((0.5 - a_h) - ((0.5 + a_h) * 2 * G / k) - ((0.25 - a_h^2) * 2 * F));

delta_R_A = (1 - (g_h * g_alpha)) * mu^2 * r_alpha^2 * freq_h^2 / freq_alpha^2;
delta_R_B = (mu * freq_h^2 * (E_R - (g_h * E_I))) + (mu * r_alpha^2 * (A_R - (g_alpha * A_I)));
delta_R_C = (A_R * E_R) - (B_R * D_R) - (A_I * E_I) + (B_I * D_I);

delta_I_A = (g_h + g_alpha) * mu^2 * r_alpha^2 * freq_h^2 / freq_alpha^2;
delta_I_B = (mu * freq_h^2 * ((g_h * E_R) + E_I)) + (mu * r_alpha^2 * (A_I + (g_alpha * A_R)));
delta_I_C = (A_I * E_R) - (B_R * D_I) + (A_R * E_I) - (B_I * D_R);

% X = freq_alpha^2 / freq^2;

% A = A_R + (i * A_I) + (1 + (i * g_h)) * mu * (freq_h / freq_alpha)^2 * X;
% B = B_R + (i * B_I);
% D = D_R + (i * D_I);
% E = E_R + (i * E_I) + (1 + (i * g_h)) * mu * r_alpha^2 * X;
% 
% Delta = det([A, B; D, E]);
