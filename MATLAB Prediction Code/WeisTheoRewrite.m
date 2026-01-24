%% Rewrite of the Weisshaar Metrix method for the theodorsen function
% Alexander Ketzle
clc, clear

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
% b = 2.25/12; % ft
% x_alpha = 0.0;
% r_alpha = 0.57785;
% freq_alpha = 2866.77844; % rad/s
% freq_h = 1853.69007; % rad/s
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

% weisshaar page 424 example
% a_h = -0.2;
% b = 3; % ft
% x_alpha = 0.1;
% r_alpha = 0.5;
% freq_h = 10;
% freq_alpha = 25;
% mu = 20;

%% Calculation

i = sqrt(-1);
kStepSize = 0.00001;
kRange = [kStepSize,6];
n = ((max(kRange) - min(kRange)) / kStepSize) + 1;
kSet = linspace(kRange(2),kRange(1),n);

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

omega(1,:) = ((A11(:) + A22(:)) - sqrt((A11(:) + A22(:)).^2 - (4 * (A11(:) .* A22(:) - A12(:) .* A21(:))))) ./ 2; 
omega(2,:) = ((A11(:) + A22(:)) + sqrt((A11(:) + A22(:)).^2 - (4 * (A11(:) .* A22(:) - A12(:) .* A21(:))))) ./ 2;

omega_r = real(omega);
omega_i = imag(omega);
freq_f = freq_alpha ./ real(sqrt(omega_r));
V_f = (freq_f .* b ./ k);
eigRatio = (omega_i ./ omega_r);
freqRatio = freq_f ./ freq_alpha;

%% Plotting

if(max(eigRatio(1,:)) > 0)
    [rowf,colf] = find(abs(eigRatio(1,:)) == min(abs(eigRatio(1,:))));
    Uf = V_f(rowf,colf);
    flutterFreq = freq_f(rowf,colf) / (2 * pi());
    [~,cold] = find(abs(eigRatio(2,:)) == min(abs(eigRatio(2,:))));
    Ud = V_f(2,cold);
    freq_div = freq_f(2,cold) / (2 * pi());
else
    [rowf,colf] = find(abs(eigRatio(2,:)) == min(abs(eigRatio(2,:))));
    Uf = V_f(rowf,colf);
    flutterFreq = freq_f(rowf,colf) / (2 * pi());
    [~,cold] = find(abs(eigRatio(1,:)) == min(abs(eigRatio(1,:))));
    Ud = V_f(1,cold);
    freq_div = freq_f(1,cold) / (2 * pi());
end

fprintf("Uf = %g ft/s\nFcr = %g Hz\nUd = %g ft/s\nFd = %g hz\n",Uf,flutterFreq,Ud,freq_div)


% U vs g
figure();
hold on;
plot(V_f(1,:),eigRatio(1,:))
plot(V_f(2,:),eigRatio(2,:))
yline(0)
ylabel("Xi/Xr");
xlabel("Velocity");
xlim([0 max(V_f(2,:) * 1.2)])



% freq ratio vs U
% figure();
% hold on;
% plot(V_f(1,:),freqRatio(1,:))
% plot(V_f(2,:),freqRatio(2,:))
% xlim([0 max(V_f(2,:) * 1.2)])
% min(freqRatio(1,:))