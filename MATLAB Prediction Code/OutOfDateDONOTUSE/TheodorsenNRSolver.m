%% Theodorsen Newton-Raphsen solver
% Alexander Ketzle
clc, clear;
[X1,k1] = TheodorsenNRSolver2(1e-6,1e2,1,0.001)

function [X,k] = TheodorsenNRSolver2(convergence, iterations, X0, k0)
    %% Given values
    freq_alpha = 4;
    freq_h = 1;
    r_alpha = sqrt(0.25);
    r_h = 0;
    b = 1;
    c = 0.6;
    rho = 0;
    m = 0;
    a = -0.4;
    x_alpha = 0.2;
    r_beta = sqrt(0.0012);
    freq_beta = 1.22;
    x_beta = 0;
    g_alpha = 0;
    g_beta = 0;
    g_h = 0;
    kappa = 0.25;
    
    %% True Constants
    p = -(1/3) * (sqrt(1 - c^2))^3;
    T_1 = -(1/3) * sqrt(1 - c^2) * (2 + c^2) + (c * acos(c));
    T_2 = c * (1 - c^2) - sqrt(1 - c^2) * (1 + c^2) * acos(c) + c*acos(c)^2;
    T_3 = -((1/8) + c^2) * acos(c)^2 + (1/4) * c * sqrt(1 - c^2) * acos(c) * (7 + (2 * c^2)) - ((1/8) * (1 - c^2) * (5 * c^2 + 4));
    T_4 = -acos(c) + c * sqrt(1- c^2);
    T_5 = -(1 - c^2) - acos(c)^2 + (2 * c * sqrt(1 - c^2) * acos(c));
    T_6 = T_2;
    T_7 = -((1/8) + c^2) * acos(c) + ((1/8) * c * sqrt(1 - c^2) * (7 + 2 * c^2));
    T_8 = -((1/3) * sqrt(1 - c^2) * (2 * c^2 + 1)) + c * acos(c);
    T_9 = (1/2) * (((1/3) * sqrt(1 - c^2)) + (a * T_4)); % WARNING: there's an inconsistency in TR496 about this one, look into it!!!
    T_10 = sqrt(1 - c^2) + acos(c);
    T_11 = acos(c)* (1 - 2 * c) + (sqrt(1 - c^2) * (2 - c));
    T_12 = sqrt(1 - c^2) * (2 + c) - (acos(c) * (2 * c + 1));
    T_13 = (1/2) * (-T_7 - ((c - a) * T_1));
    T_14 = (1/16) + (a * c / 2);

    A_alpha1 = (r_alpha^2 / kappa) + ((1/8) + a^2);
    A_alpha2 = 0.5 - a;
    A_beta1 = (r_beta^2 / kappa) - ((T_7 / pi())) + ((c - a) * ((x_beta / kappa) - (T_1 / pi())));
    A_beta2 = (1 / pi()) * (-2 * p - ((1/2) - a) * T_1);
    A_beta3 = (T_4 + T_10) / pi();
    A_h1 = (x_alpha / kappa) - a;
    B_alpha1 = A_beta1;
    B_alpha2 = (p - T_1 - (T_4 / 2)) / pi();
    B_beta1 = (r_beta^2 / kappa) - (T_3 / pi()^2);
    B_beta2 = -(T_4 * T_11) / (2 * pi()^2);
    B_beta3 = (T_5 - (T_4 * T_10)) / pi()^2;
    B_h1 = (x_beta / kappa) - (T_1 / pi());
    C_alpha1 = A_h1;
    C_alpha2 = 1;
    C_beta1 = B_h1;
    C_beta2 = - T_4 / pi();
    C_beta3 = 0;
    C_h1 = 1 + (1 / kappa);

    A_1 = det([A_alpha1, A_h1; C_alpha1, C_h1]);
    B_1 = det([A_alpha1, -(0.5 + a); C_alpha1, 1]) + ((0.5 - a) * det([-(0.5 + a), A_h1; 1, C_h1]));
    C_1 = det([(A_h1 - A_alpha2), -(0.5 + a); (C_h1 - C_alpha2), 1]);
    D_1 = -det([A_alpha2, A_h1; C_alpha2, C_h1]);

    %% Dependent functions
    sympref("FloatingPointOutput",true);
    syms k X;
    
    J_1 = besselj(1,k);
    J_0 = besselj(0,k);
    Y_0 = bessely(0,k);
    Y_1 = bessely(1,k);
    
    F = (J_1 * (J_1 + Y_0) + (Y_1 * (Y_1 - J_0))) / ((J_1 +Y_0)^2 + (Y_1 - J_0)^2);
    G = - ((Y_1 * Y_0) + (J_1 * J_0)) / ((J_1 + Y_0)^2 + (Y_1 - J_0)^2);
    
    % I and R functions
    
    R_aalpha = -A_alpha1 + ((0.25 - a^2) * (2 * G / k)) - ((0.5 + a) * (2 * F / k^2));
    R_abeta = -A_beta1 + (A_beta3 / k^2) + (0.5 + a) * ((T_11 * G / (pi() * k)) - (T_10 * F * 2 / (pi() * k^2)));
    R_ah = -A_h1 + (0.5 + a) * (2 * G / k);
    
    R_balpha = -B_alpha1 - (T_12 * (((0.5 - a) * 2 * G / k) - (2 * F / k^2)) / (2 * pi()));
    R_bbeta = -B_beta1 + (B_beta3 / k^2) - (T_12 * ((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / (pi() * k^2))) / (2 * pi()));
    R_bh = -B_h1 - (T_12 * G / k) / pi();
    
    R_calpha = -C_alpha1 - ((0.5 - a) * 2 * G / k) + (2 * F / k^2);
    R_cbeta = -C_beta1 - (T_11 * G / (pi() * k)) + (T_10 * 2 * F / (pi() * k^2));
    R_ch = -C_h1 - (2 * G / k);
    
    I_aalpha = (1/k) * (-(((1/2) + a) * 2 * G / k) - (((1/4) - a^2) * 2 * F) + A_alpha2);
    I_abeta = (-(0.5 + a) * (((T_10 * G * 2) / (pi() * k)) + (T_11 * F / pi())) + A_beta2) / k;
    I_ah = -(0.5 + a) * 2 * F;
    
    I_balpha = ((T_12 / (2 * pi())) * ((2 * G / k) + ((0.5 - a) * (2 * F))) + B_alpha2) / k;
    I_bbeta = (((T_12 / (2 * pi())) * ((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / pi()))) + B_beta2) / k;
    I_bh = (T_12 * F / pi()) / k;
    
    I_calpha = ((2 * G / k) - ((0.5 - a) * 2 * F) + C_alpha2) / k;
    I_cbeta = (((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / pi())) + C_beta2) / k;
    I_ch = (1/k) * 2 * F;
    
    % Omega functions
    omega_h = (freq_h / freq_alpha)^2 / r_alpha^2;
    omega_alpha = 1;
    
    % more coeffs for calculation
    M_1Real = A_1 + (2 * G * B_1 / k) + (2 * C_1 * F / k^2);
    Coeff_X2Real = omega_h * omega_alpha * (1 - (g_h * g_alpha));
    Coeff_XReal = (omega_h * (R_aalpha - (g_h * I_aalpha))) + (omega_alpha * (R_ch - (g_alpha * I_ch)));
    
    M_1Imag = (1/k) * (D_1 + (2 * G * C_1 / k) - (B_1 * 2 * F));
    Coeff_X2Imag = omega_h * omega_alpha * (g_h + g_alpha);
    Coeff_XImag = omega_h * ((R_aalpha * g_h) + I_aalpha) + (omega_alpha * ((R_ch * g_alpha) + I_ch));

    x_realeq = (Coeff_X2Real * X^2) + (Coeff_XReal * X) + M_1Real;
    x_imageq = (Coeff_X2Imag * X^2) + (Coeff_XImag * X) + M_1Imag;

    for i=1:iterations
        fprintf("Iteration: %d\n     X0: %f\n     k0: %f\n",i,X0,k0);
        J_1 = besselj(1,k);
        J_0 = besselj(0,k);
        Y_0 = bessely(0,k);
        Y_1 = bessely(1,k);
        
        F = (J_1 * (J_1 + Y_0) + (Y_1 * (Y_1 - J_0))) / ((J_1 +Y_0)^2 + (Y_1 - J_0)^2);
        G = - ((Y_1 * Y_0) + (J_1 * J_0)) / ((J_1 + Y_0)^2 + (Y_1 - J_0)^2);
        
        % I and R functions
        
        R_aalpha = -A_alpha1 + ((0.25 - a^2) * (2 * G / k)) - ((0.5 + a) * (2 * F / k^2));
        R_abeta = -A_beta1 + (A_beta3 / k^2) + (0.5 + a) * ((T_11 * G / (pi() * k)) - (T_10 * F * 2 / (pi() * k^2)));
        R_ah = -A_h1 + (0.5 + a) * (2 * G / k);
        
        R_balpha = -B_alpha1 - (T_12 * (((0.5 - a) * 2 * G / k) - (2 * F / k^2)) / (2 * pi()));
        R_bbeta = -B_beta1 + (B_beta3 / k^2) - (T_12 * ((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / (pi() * k^2))) / (2 * pi()));
        R_bh = -B_h1 - (T_12 * G / k) / pi();
        
        R_calpha = -C_alpha1 - ((0.5 - a) * 2 * G / k) + (2 * F / k^2);
        R_cbeta = -C_beta1 - (T_11 * G / (pi() * k)) + (T_10 * 2 * F / (pi() * k^2));
        R_ch = -C_h1 - (2 * G / k);
        
        I_aalpha = (1/k) * (-(((1/2) + a) * 2 * G / k) - (((1/4) - a^2) * 2 * F) + A_alpha2);
        I_abeta = (-(0.5 + a) * (((T_10 * G * 2) / (pi() * k)) + (T_11 * F / pi())) + A_beta2) / k;
        I_ah = -(0.5 + a) * 2 * F;
        
        I_balpha = ((T_12 / (2 * pi())) * ((2 * G / k) + ((0.5 - a) * (2 * F))) + B_alpha2) / k;
        I_bbeta = (((T_12 / (2 * pi())) * ((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / pi()))) + B_beta2) / k;
        I_bh = (T_12 * F / pi()) / k;
        
        I_calpha = ((2 * G / k) - ((0.5 - a) * 2 * F) + C_alpha2) / k;
        I_cbeta = (((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / pi())) + C_beta2) / k;
        I_ch = (1/k) * 2 * F;
        
        % Omega functions
        omega_h = (freq_h / freq_alpha)^2 / r_alpha^2;
        omega_alpha = 1;
        
        % more coeffs for calculation
        M_1Real = A_1 + (2 * G * B_1 / k) + (2 * C_1 * F / k^2);
        Coeff_X2Real = omega_h * omega_alpha * (1 - (g_h * g_alpha));
        Coeff_XReal = (omega_h * (R_aalpha - (g_h * I_aalpha))) + (omega_alpha * (R_ch - (g_alpha * I_ch)));
        
        M_1Imag = (1/k) * (D_1 + (2 * G * C_1 / k) - (B_1 * 2 * F));
        Coeff_X2Imag = omega_h * omega_alpha * (g_h + g_alpha);
        Coeff_XImag = omega_h * ((R_aalpha * g_h) + I_aalpha) + (omega_alpha * ((R_ch * g_alpha) + I_ch));
    
        x_realeq = (Coeff_X2Real * X^2) + (Coeff_XReal * X) + M_1Real == 0;
        x_imageq = (Coeff_X2Imag * X^2) + (Coeff_XImag * X) + M_1Imag == 0;
        
        [Xi, ki] = vpasolve([x_realeq, x_imageq],[X, k],[X0, k0]);
        Xdiff = Xi - X0;
        kdiff = ki - k0;

        if(Xdiff < convergence && kdiff < convergence)
            fprintf("Converged in %d iterations\n",i);
            break
        elseif(abs(Xdiff) > 10 && abs(kdiff) > 10)
            fprintf("Both solutions diverging, attempting correction\n");
            fprintf("     X solution:      %f\n",Xi);
            fprintf("     X solution diff: %f\n",Xdiff);
            fprintf("     k solution:      %f\n",ki);
            fprintf("     k solution diff: %f\n",kdiff);
            k0 = k0 / 2;
            X0 = X0 / 2;
        elseif(abs(kdiff) > 10)
            fprintf("k Solution diverging, attempting correction\n");
            fprintf("     k solution:      %f\n",ki);
            fprintf("     k solution diff: %f\n",kdiff);
            X0 = X0 / 2;
            k0 = k0 / 2;
        elseif(abs(Xdiff) > 10)
            fprintf("X Solution diverging, attempting correction\n");
            fprintf("     X solution:      %f\n",Xi);
            fprintf("     X solution diff: %f\n",Xdiff);
            X0 = X0 / 2;
        else
            X0 = Xi;
            k0 = ki;
        end
    end
    fprintf("Calculated Solution:\n     X = %f\n     k = %f\n",Xi,ki);
    X = Xi;
    k = ki;
end