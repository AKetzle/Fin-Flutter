%% Weisshaar Modified Theodorsen Method
% Alexander Ketzle
% Long Suffering
clc, clear, close all;

% TR685 parameters - should converge near 831.6 ft/s and k = 0.406
% a = -0.4;
% b = 6; % ft
% x_thetabar = 0.2;
% r_thetabar = sqrt(0.25);
% freq_theta = 90; % rad/s
% freq_h = freq_theta / 4; % rad/s
% mu = 4;

% Weisshaar parameters - should converge near 166 ft/s and 216 ft/s
% depending on flutter vs divergence - depends on eigenvalue used
% a = -0.2;
% b = 3; % ft
% x_thetabar = 0.1;
% r_thetabar = 0.5;
% freq_theta = 25; % rad/s
% freq_h = 10; % rad/s
% mu = 20;

% Cippola N5800 Parameters - This is one of the examples shipped with
%FinSim
% a = 0.0;
% b = 3.5625 / 12; % ft
% x_thetabar = 0.316;
% r_thetabar = 0.57757;
% freq_theta = 2593.373; % rad/s
% freq_h = 2458.08525; % rad/s
% mu = 77.11441;

% Y.C. Fung Parameters - Page 236
% a = 0.0;
% b = 30; % ft
% x_thetabar = 0.0;
% r_thetabar = sqrt(0.6222);
% freq_theta = sqrt(2.41);
% freq_h = sqrt(0.775);
% mu = 40;

% Y.C. Fung parameters - Page 219
% a =  -0.15;
% b = 2;
% x_thetabar = 0.25;
% r_thetabar = sqrt(0.388);
% freq_theta = 48;
% freq_h = (55.9 / 64.1) * 48;
% mu = 76;
<<<<<<< Updated upstream
<<<<<<< Updated upstream
=======
=======
>>>>>>> Stashed changes

% Doppler's Fin
a =  0.0;
b = 2.25/12;
x_thetabar = 0.0;
r_thetabar = 0.57785;
freq_theta = 2866.77844;
freq_h = 1853.69007;
mu = 110.28942;
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

i = sqrt(-1);
% I'm keeping these equations here so I don't lose them
% syms h_bar theta_bar reduced_vel flutter_freq k
% 
% J_1 = besselj(1,k);
% Y_1 = bessely(1,k);
% Y_0 = bessely(0,k);
% J_0 = besselj(0,k);
% 
% C = (J_1 - (i * Y_1)) / ((J_1 + Y_0) + (i * (J_0 - Y_1)));
% 
% L_h = 1 - (i * 2 * C / k);
% L_alpha = 0.5 - (i * (1 + (2 * C) / k)) - (2 * C / k^2);
% M_alpha = (3/8) - (i / k);
% 
% fluttermatrix = [(freq_theta^2 / freq_h^2) * (1 + (L_h / mu)), (freq_theta^2 / freq_h^2) * (x_thetabar + (mu^-1 * (L_alpha - ((0.5 + a) * L_h))));
%     ((x_thetabar / r_thetabar^2) + ((1 / (mu * r_thetabar^2)) * (0.5 - (0.5 + a) * L_h))), 1 + ((1 / (mu * r_thetabar^2)) * (M_alpha - ((L_alpha + 0.5) * (0.5 + a)) + (L_h * (0.5 + a)^2)))];
% omega = eig(fluttermatrix);

%% New Solver Attempt

% set up conditions

velStepSize = 1; % ft/s per step
<<<<<<< Updated upstream
<<<<<<< Updated upstream
vel_range = [50,2000]; % ft/s, range of values to test
=======
vel_range = [50,3000]; % ft/s, range of values to test
>>>>>>> Stashed changes
=======
vel_range = [50,3000]; % ft/s, range of values to test
>>>>>>> Stashed changes
n = ((vel_range(2) - vel_range(1)) / velStepSize) + 1;
testVels = linspace(vel_range(1), vel_range(2), n);
solutionMatrix1 = zeros([6,size(testVels,2)]); % each column corresponds to a test velocity


initial_k = 0.226;
stepLimiter = 1.0; % scaling factor to reduce how much value difference comes into play in solver
iterations = 50;
convergence = 0.0001;
whicheigenvalue = 1;
% initialize the function

parfor velStep = 1:n
    V = testVels(velStep);
    k = initial_k;
    for iter = 1:iterations

        J_0 = besselj(0,k);
        J_1 = besselj(1,k);
        Y_0 = bessely(0,k);
        Y_1 = bessely(1,k);
        
        C = ((i * Y_1) - J_1) / (-(J_1 + Y_0) + (i * (Y_1 - J_0)));
        
        L_h = 1 - (i * 2 * C / k);
        L_alpha = 0.5 - (i * (1 + (2 * C)) / k) - (2 * C / k^2);
        M_alpha = (3/8) - (i / k);
        
        A11 = (freq_theta^2 / freq_h^2) * (1 + (L_h / mu));
        A12 = (freq_theta^2 / freq_h^2) * (x_thetabar + (mu^-1 * (L_alpha - ((0.5 + a) * L_h))));
        A21 = ((x_thetabar / r_thetabar^2) + ((1 / (mu * r_thetabar^2)) * (0.5 - (0.5 + a) * L_h)));
        A22 = 1 + ((1 / (mu * r_thetabar^2)) * (M_alpha - ((L_alpha + 0.5) * (0.5 + a)) + (L_h * (0.5 + a)^2)));

        fluttermatrix = [A11, A12; A21, A22];
        omega = eig(fluttermatrix);

        eigen = subs(omega); % calculate function
        freq_f = freq_theta / sqrt(real(eigen(whicheigenvalue))); % extract necessary component from function
        k1 = double(freq_f * b / V); % calculate new value
        w = k1 - k; % calculate residual
        if abs(w) > convergence && iter < iterations
            k = k + (stepLimiter * w); % produce next value to test
        else
            V_f = freq_f * b / k;
            fprintf("Calculation Complete in %d iterations\nFor V = %g ft/s:\nk = %g\nV_f = %g ft/s\n",iter,V,k,V_f);
            solutionMatrix1(:,velStep) = [V,V_f,k,eigen(1),eigen(2),freq_f].';
            break;
        end
    end
end



solutionMatrix2 = zeros([6,size(testVels,2)]); % each column corresponds to a test velocity
whicheigenvalue = 2;
% initialize the function

parfor velStep = 1:n
    V = testVels(velStep);
    k = initial_k;
    for iter = 1:iterations

        J_0 = besselj(0,k);
        J_1 = besselj(1,k);
        Y_0 = bessely(0,k);
        Y_1 = bessely(1,k);
        
        C = ((i * Y_1) - J_1) / (-(J_1 + Y_0) + (i * (Y_1 - J_0)));
        
        L_h = 1 - (i * 2 * C / k);
        L_alpha = 0.5 - (i * (1 + (2 * C)) / k) - (2 * C / k^2);
        M_alpha = (3/8) - (i / k);
        
        A11 = (freq_theta^2 / freq_h^2) * (1 + (L_h / mu));
        A12 = (freq_theta^2 / freq_h^2) * (x_thetabar + (mu^-1 * (L_alpha - ((0.5 + a) * L_h))));
        A21 = ((x_thetabar / r_thetabar^2) + ((1 / (mu * r_thetabar^2)) * (0.5 - (0.5 + a) * L_h)));
        A22 = 1 + ((1 / (mu * r_thetabar^2)) * (M_alpha - ((L_alpha + 0.5) * (0.5 + a)) + (L_h * (0.5 + a)^2)));

        fluttermatrix = [A11, A12; A21, A22];
        omega = eig(fluttermatrix);

        eigen = subs(omega); % calculate function
        freq_f = freq_theta / sqrt(real(eigen(whicheigenvalue))); % extract necessary component from function
        k1 = double(freq_f * b / V); % calculate new value
        w = k1 - k; % calculate residual
        if abs(w) > convergence && iter < iterations
            k = k + (stepLimiter * w); % produce next value to test
        else
            V_f = freq_f * b / k;
            fprintf("Calculation Complete in %d iterations\nFor V = %g ft/s:\nk = %g\nV_f = %g ft/s\n",iter,V,k,V_f);
            solutionMatrix2(:,velStep) = [V,V_f,k,eigen(1),eigen(2),freq_f].';
            break;
        end
    end
end

%% - plotting

freqRatio1 = imag(solutionMatrix1(4,:)) ./ real(solutionMatrix1(4,:));
freqRatio2 = imag(solutionMatrix2(4,:)) ./ real(solutionMatrix2(4,:));
freqRatio3 = imag(solutionMatrix1(5,:)) ./ real(solutionMatrix2(5,:));
freqRatio4 = imag(solutionMatrix2(5,:)) ./ real(solutionMatrix2(5,:));
figure;
hold on;
grid on;
plot(testVels,freqRatio1)
plot(testVels,freqRatio2)
plot(testVels,freqRatio3)
plot(testVels,freqRatio4)
legend("First Freq Ratio","Second Freq Ratio","Third Freq Ratio","Fourth Freq Ratio","AutoUpdate","off")
yline(0)
ylabel("Xi/Xr");
xlabel("velocity");
%axis([0 inf -1 0.2])
%xticks(0:20:400);
hold off;

[val1, idx1] = find(abs(freqRatio1(:,:)) == min(abs(freqRatio1)));
solutionMatrix1(:,idx1)
[val2, idx2] = find(abs(freqRatio2(:,:)) == min(abs(freqRatio2)));
solutionMatrix2(:,idx2)
[val3, idx3] = find(abs(freqRatio3(:,:)) == min(abs(freqRatio3)));
solutionMatrix1(:,idx3)
[val4, idx4] = find(abs(freqRatio4(:,:)) == min(abs(freqRatio4)));
solutionMatrix2(:,idx4)

thetafreqratio1 = solutionMatrix1(6,:) ./ freq_theta;
thetafreqratio2 = solutionMatrix2(6,:) ./ freq_theta;
figure;
hold on;
grid on;
plot(testVels,thetafreqratio1)
plot(testVels,thetafreqratio2)
legend("First Freq Ratio","Second Freq Ratio","AutoUpdate","off")
ylabel("freq_f/freq_theta");
xlabel("velocity");
%axis([0 inf 0 1.2])
hold off;

sqrtXr1 = freq_theta ./ solutionMatrix1(6,:);
sqrtXr2 = freq_theta ./ solutionMatrix2(6,:);
sqrtXi1 = sqrt(-i .* imag(solutionMatrix1(4,:)));
sqrtXi2 = sqrt(-i .* imag(solutionMatrix2(4,:)));
invK1 = 1 ./ solutionMatrix1(3,:);
invK2 = 1 ./ solutionMatrix2(3,:);
figure;
hold on;
grid on;
plot(invK1,sqrtXr1)
plot(invK1,1 ./ thetafreqratio2)
plot(invK1,sqrtXi1)
plot(invK1,sqrtXi2)
xlabel("1/k");
ylabel("sqrt(X)");
legend("sqrt(Xr1)","sqrt(Xr2)","sqrt(Xi1)","sqrt(Xi2)");
axis([0 6 0.9 2.6])
xticks(0:0.6:6)
yticks(0.9:0.17:2.6)
hold off;