%% Y.C. Fung Ch. 6.9 Flutter Calculation
% based on Theodorsen Method
% Alexander Ketzle
% starts around page 213
clc, clear;
i = sqrt(-1);

% TR685 parameters - should converge near 831.6 ft/s and k = 0.406
% a_h = -0.4;
% b = 6; % ft
% x_alpha = 0.2;
% r_alpha = sqrt(0.25);
% freq_alpha = 90; % rad/s
% freq_h = freq_alpha / 4; % rad/s
% mu = 4;
% g_h = 0;
% g_alpha = 0;


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

% page 236 parameters
% mu = 40;
% a_h = 0.0;
% x_alpha = 0.0;
% r_alpha = sqrt(0.6222);
% b = 30; % ft
% freq_alpha = sqrt(2.41); % rad/s
% freq_h = sqrt(0.775); % rad/s
% g_alpha = 0;
% g_h = 0;

% Cippola N5800 Parameters - This is one of the examples shipped with
%FinSim
a_h = 0.0;
b = 3.5625 / 12; % ft
x_alpha = 0.316;
r_alpha = 0.57757;
freq_alpha = 2593.373; % rad/s
freq_h = 2458.08525; % rad/s
mu = 77.11441;
g_alpha = 0.005;
g_h = 0.005;


invkstepsize = 0.00001;
invkrange = [invkstepsize,6.4];
n = uint32(((invkrange(2) - invkrange(1)) / invkstepsize) + 1);
invk_set = linspace(invkrange(1),invkrange(2),n);
solutionmatrix = zeros([13,size(invk_set,2)]);
imagMatrix1 = zeros([3,size(invk_set,2)]);

parfor invkstep = 1:n
    k_inv = invk_set(invkstep);
    k = 1 / k_inv;

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
    delta_R_B = (mu * freq_h^2 / freq_alpha^2 * (E_R - (g_h * E_I))) + (mu * r_alpha^2 * (A_R - (g_alpha * A_I)));
    delta_R_C = (A_R * E_R) - (B_R * D_R) - (A_I * E_I) + (B_I * D_I);
    
    delta_I_A = (g_h + g_alpha) * mu^2 * r_alpha^2 * freq_h^2 / freq_alpha^2;
    delta_I_B = (mu * freq_h^2 / freq_alpha^2 * ((g_h * E_R) + E_I)) + (mu * r_alpha^2 * (A_I + (g_alpha * A_R)));
    delta_I_C = (A_I * E_R) - (B_R * D_I) + (A_R * E_I) - (B_I * D_R);

    X_R1 = (-delta_R_B + sqrt(delta_R_B^2 - (4 * delta_R_A * delta_R_C))) / (2 * delta_R_A);
    X_R2 = (-delta_R_B - sqrt(delta_R_B^2 - (4 * delta_R_A * delta_R_C))) / (2 * delta_R_A);
    X_I1 = -delta_I_C / delta_I_B;
    X_I2 = -delta_I_C / delta_I_B;
    rt_X_R1 = sqrt(X_R1);
    rt_X_R2 = sqrt(X_R2);
    rt_X_I1 = sqrt(X_I1);
    rt_X_I2 = sqrt(X_I2);
    imagMatrix1(:,invkstep) = [delta_I_A; delta_I_B; delta_I_C];
    realMatrix1(:,invkstep) = [delta_R_A; delta_R_B; delta_R_C];
    solutionmatrix(:,invkstep) = [k, k_inv, X_R1, X_R2, X_I1, X_I2, rt_X_R1, rt_X_R2, rt_X_I1, rt_X_I2,delta_I_A,delta_I_B,delta_I_C].';
end

%% Plotting

XRatio1 = abs(1 - abs(solutionmatrix(7,:) ./ solutionmatrix(9,:)));
XRatio2 = abs(1 - abs(solutionmatrix(8,:) ./ solutionmatrix(9,:)));


figure();
hold on;
invk = solutionmatrix(2,:);
rt_X = solutionmatrix(7:10,:);
plot(invk,rt_X)
xlabel("1/k");
ylabel("sqrt(X)");


[~,idx1] = find(XRatio2 == min(XRatio2,[],"omitnan"))
flutterPoint = solutionmatrix(:,idx1)
Uf = freq_alpha * b * flutterPoint(2) / flutterPoint(8)
flutterFreq = freq_alpha / flutterPoint(8)
inv_kf = flutterPoint(2)
kf = flutterPoint(1)
flutter_rtX = flutterPoint(8)

% X = freq_alpha^2 / freq^2;
% A = A_R + (i * A_I) + (1 + (i * g_h)) * mu * (freq_h / freq_alpha)^2 * X;
% B = B_R + (i * B_I);
% D = D_R + (i * D_I);
% E = E_R + (i * E_I) + (1 + (i * g_h)) * mu * r_alpha^2 * X;
% Delta = det([A, B; D, E]);