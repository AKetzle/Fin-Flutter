%% Weisshaar Modified Theodorsen Method
% Alexander Ketzle
% Long Suffering
clc, clear, close all;

% TR685 parameters - should converge near 831.6 ft/s and k = 0.406
a = -0.4;
b = 6; % ft
x_thetabar = 0.2;
r_thetabar = sqrt(0.25);
freq_theta = 90; % rad/s
freq_h = freq_theta / 4; % rad/s
mu = 0.25;

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
% FinSim
% a = 0.0;
% b = 2 * 3.5625 / 12; % ft
% x_thetabar = 0.1;
% r_thetabar = 0.57757;
% freq_theta = 2593.373; % rad/s
% freq_h = 2458.08525; % rad/s
% mu = 77.11441;
% initial_k = 0.428


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

velStepSize = 0.1; % ft/s per step
vel_range = [200,1300]; % ft/s, range of values to test
n = ((vel_range(2) - vel_range(1)) / velStepSize) + 1;
testVels = linspace(vel_range(1), vel_range(2), n);
solutionMatrix = zeros([4,size(testVels,2)]); % each column corresponds to a test velocity


initial_k = 0.406;
stepLimiter = 1.0; % scaling factor to reduce how much value difference comes into play in solver
iterations = 150;
convergence = 0.0001;
% initialize the function

parfor velStep = 1:n
    V = testVels(velStep);
    k = initial_k;
    for iter = 1:iterations

        J_1 = besselj(1,k);
        Y_1 = bessely(1,k);
        Y_0 = bessely(0,k);
        J_0 = besselj(0,k);
        
        C = (J_1 - (i * Y_1)) / ((J_1 + Y_0) + (i * (J_0 - Y_1)));
        
        L_h = 1 - (i * 2 * C / k);
        L_alpha = 0.5 - (i * (1 + (2 * C) / k)) - (2 * C / k^2);
        M_alpha = (3/8) - (i / k);
        

        % fluttermatrix = [(freq_theta^2 / freq_h^2) * (1 + (L_h / mu)), (freq_theta^2 / freq_h^2) * (x_thetabar + (mu^-1 * (L_alpha - ((0.5 + a) * L_h))));
        %     ((x_thetabar / r_thetabar^2) + ((1 / (mu * r_thetabar^2)) * (0.5 - (0.5 + a) * L_h))), 1 + ((1 / (mu * r_thetabar^2)) * (M_alpha - ((L_alpha + 0.5) * (0.5 + a)) + (L_h * (0.5 + a)^2)))];
        A11 = (freq_theta^2 / freq_h^2) * (1 + (L_h / mu));
        A12 = (freq_theta^2 / freq_h^2) * (x_thetabar + (mu^-1 * (L_alpha - ((0.5 + a) * L_h))));
        A21 = ((x_thetabar / r_thetabar^2) + ((1 / (mu * r_thetabar^2)) * (0.5 - (0.5 + a) * L_h)));
        A22 = 1 + ((1 / (mu * r_thetabar^2)) * (M_alpha - ((L_alpha + 0.5) * (0.5 + a)) + (L_h * (0.5 + a)^2)));

        fluttermatrix = [A11, A12; A21, A22];
        omega = eig(fluttermatrix);

        eigen = subs(omega); % calculate function
        freq_f = freq_theta / sqrt(real(eigen(1))); % extract necessary component from function
        k1 = double(freq_f * b / V); % calculate new value
        w = k1 - k; % calculate residual
        if abs(w) > convergence && iter < iterations
            k = k + (stepLimiter * w); % produce next value to test
        else
            V_f = freq_f * b / k;
            fprintf("Calculation Complete in %d iterations\nFor V = %g ft/s:\nk = %g\nV_f = %g ft/s\n",iter,V,k,V_f);
            solutionMatrix(:,velStep) = [k,V_f,eigen(1),eigen(2)].';
            break;
        end
    end
end

freqRatio1 = imag(solutionMatrix(3,:)) ./ real(solutionMatrix(3,:));
freqRatio2 = imag(solutionMatrix(4,:)) ./ real(solutionMatrix(4,:));
hold on;
plot(testVels,freqRatio1)
plot(testVels,freqRatio2)
legend("First Freq Ratio","Second Freq Ratio","AutoUpdate","off")
yline(0)
[val, idx] = find(abs(imag(solutionMatrix(4,:)) ./ real(solutionMatrix(4,:))) == min(abs(freqRatio2)));
solutionMatrix(:,idx)
