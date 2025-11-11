%% Weisshaar Modified Theodorsen Method
% Alexander Ketzle
% Long Suffering
clc, clear, close all;

a = -0.2;
b = 3; % ft
x_thetabar = 0.1;
r_thetabar = 0.5;
freq_h = 10; % rad/s
freq_theta = 25; % rad/s
mu = 20;
i = sqrt(-1);
initial_guess = 1;
iterations = 10;
convergence = 1e-6;


syms h_bar theta_bar reduced_vel flutter_freq k
%k = initial_guess;
%for iter = 1:iterations
%k = reduced_freq * b / reduced_vel;    
%hold on;
%for k = 0.001:0.001:3

%end

J_1 = besselj(1,k);
Y_1 = bessely(1,k);
Y_0 = bessely(0,k);
J_0 = besselj(0,k);

C = (J_1 - (i * Y_1)) / ((J_1 + Y_0) + (i * (J_0 - Y_1)));

L_h = 1 - (i * 2 * C / k);
L_alpha = 0.5 - (i * (1 + (2 * C) / k)) - (2 * C / k^2);
M_alpha = (3/8) - (i / k);

fluttermatrix = [(freq_theta^2 / freq_h^2) * (1 + (L_h / mu)), (freq_theta^2 / freq_h^2) * (x_thetabar + (mu^-1 * (L_alpha - ((0.5 + a) * L_h))));
    ((x_thetabar / r_thetabar^2) + ((1 / (mu * r_thetabar^2)) * (0.5 - (0.5 + a) * L_h))), 1 + ((1 / (mu * r_thetabar^2)) * (M_alpha - ((L_alpha + 0.5) * (0.5 + a)) + (L_h * (0.5 + a)^2)))];
fprime = diff(fluttermatrix,k);
omega = eig(fluttermatrix);
%% N-R
k = 4.0779;
dk = 0.001;
%for iter = 1:iterations
    %iter;
    f = subs(imag(omega(1)))
    k = k + dk
    df = (imag(subs(omega(1))) - f) / dk
    k = k - dk
    term2 = (f / df)
    k1 = k - term2
    fprintf("===========\n k = %g\n i = %g\n f = %g\ndf = %g\nt2 = %g\nk1 = %g\n",k,0,f,df,term2,k1);
%     if(k1 == k)
%         break
%     else
%         k = k1;
%     end
% end
subs(omega)

freqratio = subs(1 ./ sqrt(real(omega)))
eigratio = subs(imag(omega(:,1)) ./ real(omega))
V_f = subs(freq_theta * b ./ (k * sqrt(real(omega)))) % for the givens, flutter should be 166.1 ft/s and divergence 216 ft/s
%scatter(V_f,eigratio)
%scatter(V_f,freqratio)
%yline(0)
%xlim([0,250])
%ylim([-1,0.3])



%real(omega)
%imag(omega(:,1))
%imag(omega(:,1)) / real(omega)
%flutterfreqratio = 1 ./ sqrt(real(omega)) % need elementwise math to get both flutter and divergence at once
%end