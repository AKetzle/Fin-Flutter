%% Fung Ch6 Rewrite
clc, clear;
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
% a_h = 0.0;
% b = 3.5625 / 12; % ft
% x_alpha = 0.316;
% r_alpha = 0.57757;
% freq_alpha = 2593.373; % rad/s
% freq_h = 2458.08525; % rad/s
% mu = 77.11441;
% g_alpha = 0.005;
% g_h = 0.005;


invkstepsize = 0.001;
invkrange = [invkstepsize,8];
n = uint32(((invkrange(2) - invkrange(1)) / invkstepsize) + 1);
invk_set = linspace(invkrange(1),invkrange(2),n);
k_set = 1 ./ invk_set;


k = k_set;
k_inv = invk_set;

Ch_k = besselh(1,2,k) ./ (besselh(1,2,k) + (i * besselh(0,2,k))); % Theodorsen Function; Less lines to compute than using the other bessel functions
F = real(Ch_k);
G = imag(Ch_k);

A_R = -(mu + 1) - (2 .* G ./ k);
A_I = 2 .* F ./ k;
B_R = -((mu .* x_alpha) - a_h) + (2 .* F ./ k.^2) - ((0.5 - a_h) .* 2 .* G ./ k);
B_I = (1 ./ k) .* (1 + (2 .* G ./ k) + ((0.5 - a_h) .* 2 .* F));
D_R = -((mu .* x_alpha) - a_h) + ((0.5 + a_h) .* 2 .* G ./ k);
D_I = -(0.5 + a_h) .* 2 .* F ./ k;
E_R = -((mu .* r_alpha.^2) + a_h.^2 + 0.125) + ((0.25 - a_h.^2) .* 2 .* G ./ k) - ((0.5 + a_h) .* 2 .* F ./ k.^2);
E_I = (1 ./ k) .* ((0.5 - a_h) - ((0.5 + a_h) * 2 .* G ./ k) - ((0.25 - a_h.^2) * 2 .* F));

delta_R_A = (1 - (g_h * g_alpha)) * mu^2 * r_alpha.^2 * freq_h.^2 / freq_alpha.^2;
delta_R_B = (mu * freq_h.^2 / freq_alpha.^2 .* (E_R - (g_h .* E_I))) + (mu * r_alpha.^2 .* (A_R - (g_alpha .* A_I)));
delta_R_C = (A_R .* E_R) - (B_R .* D_R) - (A_I .* E_I) + (B_I .* D_I);

delta_I_A = (g_h + g_alpha) .* mu^2 .* r_alpha.^2 .* freq_h.^2 ./ freq_alpha.^2;
delta_I_B = (mu .* freq_h.^2 ./ freq_alpha.^2 .* ((g_h .* E_R) + E_I)) + (mu .* r_alpha.^2 .* (A_I + (g_alpha .* A_R)));
delta_I_C = (A_I .* E_R) - (B_R .* D_I) + (A_R .* E_I) - (B_I .* D_R);

X_R1 = (-delta_R_B + sqrt(delta_R_B.^2 - (4 .* delta_R_A .* delta_R_C))) ./ (2 .* delta_R_A);
X_R2 = (-delta_R_B - sqrt(delta_R_B.^2 - (4 .* delta_R_A .* delta_R_C))) ./ (2 .* delta_R_A);
%X_I1 = -delta_I_C ./ delta_I_B;
X_I1 = (-delta_I_B + sqrt(delta_I_B.^2 - (4 .* delta_I_A .* delta_I_C))) ./ (2 .* delta_I_A);
X_I2 = (-delta_I_B - sqrt(delta_I_B.^2 - (4 .* delta_I_A .* delta_I_C))) ./ (2 .* delta_I_A);
%X_I2 = -delta_I_C ./ delta_I_B;
rt_X_R1 = sqrt(X_R1);
rt_X_R2 = sqrt(X_R2);
rt_X_I1 = sqrt(X_I1);
rt_X_I2 = sqrt(X_I2);

solutionmatrix = [k; k_inv; X_R1; X_R2; X_I1; rt_X_R1; rt_X_R2; rt_X_I1;];
imagMatrix2 = [zeros(1,size(delta_I_B,2)); delta_I_B; delta_I_C];
realMatrix2 = [zeros(1,size(delta_I_B,2)) + delta_R_A; delta_R_B; delta_R_C];

XRatio1 = abs(1 - abs(solutionmatrix(6,:) ./ solutionmatrix(8,:)));
XRatio2 = abs(1 - abs(solutionmatrix(7,:) ./ solutionmatrix(8,:)));


fig = figure();
hold on;
invk = solutionmatrix(2,:);
rt_X = solutionmatrix(6:8,:);
plot(invk,rt_X)
xlabel("1/k");
ylabel("sqrt(X)");
%axis([-inf inf 0.6 1.3])


[~,idx1] = find(XRatio2 == min(XRatio2,[],"omitnan"));
flutterPoint = solutionmatrix(:,idx1);
Uf = freq_alpha * b * flutterPoint(2) / flutterPoint(7);
flutterFreq = freq_alpha / flutterPoint(7) / (2 * pi());
inv_kf = flutterPoint(2);
kf = flutterPoint(1);
flutter_rtX = flutterPoint(7);

fprintf("Uf = %g ft/s\nFcr = %g Hz\nkf = %g\nsqrt(Xf) = %g\n1/kf = %g\n",Uf,flutterFreq,kf,flutter_rtX,inv_kf)