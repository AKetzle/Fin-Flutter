%% Comparison of TN4197 vs. 1/k and U vs. g methods
% Alexander Ketzle
clc, clear;

%% Fin Parameters

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
h = 4/12; % ft
t = 0.195/12; % ft
C_r = 12/12; % ft
C_t = 2.25/12; % ft
Theta_LE = 22.306; % degrees
p_0 = 14.696; % psi
p = p_0; % psf
a = 1116.4; % ft/s
G_E = 3759398.37; % psi
gamma = 1.4;

invkstepsize = 0.00001;
invkMax = 14;
kStepSize = 0.00001;
kMax = 10;

%% Flutter Speed Calculation

V_f_TN4197 = TN4197(G_E,h,t,C_r,C_t,Theta_LE,p,p_0,a);
M_TN = V_f_TN4197 / a;
M_TNc = real(sqrt(M_TN^2 * (sqrt(1 - (M_TN^4 / 4)) - (M_TN^2 / 2)))); % only keeping real portion
V_f_TN4197_corr = M_TNc * a;
V_f_invk_raw = TR496TR685(freq_alpha,freq_h,a_h,x_alpha,r_alpha,b,mu,invkstepsize,invkMax,g_h,g_alpha);
M_i1 = V_f_invk_raw / a;
M_c1 = real(sqrt(M_i1^2 * (sqrt(1 - (M_i1^4 / 4)) - (M_i1^2 / 2)))); % only keeping real portion
V_f_invk_corr = a * M_c1; % may not be correct
V_f_ug_raw = UvsgNoDamp(freq_alpha, freq_h, a_h, x_alpha, r_alpha, b, mu, kStepSize,kMax);
M_i2 = V_f_ug_raw / a;
M_c2 = real(sqrt(M_i2^2 * (sqrt(1 - (M_i2^4 / 4)) - (M_i2^2 / 2)))); % only keeping real portion
V_f_ug_corr = a * M_c2; % may not be correct
V_f_bennett_raw = BennetFlutter(C_r,C_t,Theta_LE,h,t,a,G_E,p,p_0,gamma);
M_bi = V_f_bennett_raw / a;
M_bc = real(sqrt(M_bi^2 * (sqrt(1 - (M_bi^4 / 4)) - (M_bi^2 / 2)))); % only keeping real portion
V_f_bennett_corr = M_bc * a;

% Simmons application
V_f_simmons_TN = SimmonsFlutter(x_alpha, r_alpha, freq_h, freq_alpha, b, mu, M_TN);
V_f_simmons_TN_corr = SimmonsFlutter(x_alpha, r_alpha, freq_h, freq_alpha, b, mu, M_TNc);
V_f_simmons_invk_raw = SimmonsFlutter(x_alpha, r_alpha, freq_h, freq_alpha, b, mu, M_i1);
V_f_simmons_invk_corr = SimmonsFlutter(x_alpha, r_alpha, freq_h, freq_alpha, b, mu, M_c1);
V_f_simmons_ug_raw = SimmonsFlutter(x_alpha, r_alpha, freq_h, freq_alpha, b, mu, M_i2);
V_f_simmons_ug_corr = SimmonsFlutter(x_alpha, r_alpha, freq_h, freq_alpha, b, mu, M_c2);
V_f_simmons_bennett_raw = SimmonsFlutter(x_alpha, r_alpha, freq_h, freq_alpha, b, mu, M_bi);
V_f_simmons_bennett_corr = SimmonsFlutter(x_alpha, r_alpha, freq_h, freq_alpha, b, mu, M_bc);

%% Difference Calculation

% Direct Diff
TN4197_corr_diff = (V_f_TN4197_corr - V_f_TN4197) * 100 / V_f_TN4197;
simm_TN4197_diff = (V_f_simmons_TN_corr - V_f_TN4197) * 100 / V_f_TN4197;
simm_TN4197_corr_diff = (V_f_simmons_TN_corr - V_f_TN4197) * 100 / V_f_TN4197;
invk_raw_diff = (V_f_invk_raw - V_f_TN4197) * 100 / V_f_TN4197;
invk_corr_diff = (V_f_invk_corr - V_f_TN4197) * 100 / V_f_TN4197;
simm_invk_raw_diff = (V_f_simmons_invk_raw - V_f_TN4197) * 100 / V_f_TN4197;
simm_invk_corr_diff = (V_f_simmons_invk_corr - V_f_TN4197) * 100 / V_f_TN4197;
ug_raw_diff = (V_f_ug_raw - V_f_TN4197) * 100 / V_f_TN4197;
ug_corr_diff = (V_f_ug_corr - V_f_TN4197) * 100 / V_f_TN4197;
simm_ug_raw_diff = (V_f_simmons_ug_raw - V_f_TN4197) * 100 / V_f_TN4197;
simm_ug_corr_diff = (V_f_simmons_ug_corr - V_f_TN4197) * 100 / V_f_TN4197;
bennett_raw_diff = (V_f_bennett_raw - V_f_TN4197) * 100 / V_f_TN4197;
bennett_corr_diff = (V_f_bennett_corr - V_f_TN4197) * 100 / V_f_TN4197;
simm_bennett_raw_diff = (V_f_simmons_bennett_raw - V_f_TN4197) * 100 / V_f_TN4197;
simm_bennett_corr_diff = (V_f_simmons_bennett_corr - V_f_TN4197) * 100 / V_f_TN4197;

%Simmons Diff
diff2_simm_TN4197_corr = (V_f_TN4197_corr - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_TN4197_simm_corr = (V_f_simmons_TN_corr - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_invk_raw = (V_f_invk_raw - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_invk_corr = (V_f_invk_corr - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_invk_simm_raw = (V_f_simmons_invk_raw - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_invk_simm_corr = (V_f_simmons_invk_corr - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_ug_raw = (V_f_ug_raw - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_ug_corr = (V_f_ug_corr - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_ug_simm_raw = (V_f_simmons_ug_raw - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_ug_simm_corr = (V_f_simmons_ug_corr - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_bennett_raw = (V_f_bennett_raw - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_bennett_corr = (V_f_bennett_corr - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_bennett_simm_raw = (V_f_simmons_bennett_raw - V_f_simmons_TN) * 100 / V_f_simmons_TN;
diff2_simm_bennett_simm_corr = (V_f_simmons_bennett_corr - V_f_simmons_TN) * 100 / V_f_simmons_TN;

%% Results Display

fprintf("===TN4197===\n" + ...
    "Direct:              %.1f ft/s\n" + ...
    "Direct Simmons:      %.1f ft/s\n",V_f_TN4197,V_f_simmons_TN);
fprintf("" + ...
    "Prandtl-Glauert:     %.1f ft/s\n" + ...
    "P-G Simmons:         %.1f ft/s\n",V_f_TN4197_corr,V_f_simmons_TN_corr);
fprintf("===1/k===\n" + ...
    "Direct:              %.1f ft/s\n" + ...
    "Direct Simmons:      %.1f ft/s\n" + ...
    "Prandtl-Glauert:     %.1f ft/s\n" + ...
    "P-G Simmons:         %.1f ft/s\n",V_f_ug_raw,V_f_simmons_ug_raw,V_f_ug_corr,V_f_simmons_ug_corr);
fprintf("===U vs g===\n");
fprintf( ...
    "Direct:              %.1f ft/s\n" + ...
    "Direct Simmons:      %.1f ft/s\n",V_f_ug_raw,V_f_simmons_ug_raw);
fprintf( ...
    "Prandtl-Glauert:     %.1f ft/s\n" + ...
    "P-G Simmons:         %.1f ft/s\n",V_f_ug_corr,V_f_simmons_ug_corr);
fprintf("===Bennett===\n");
fprintf( ...
    "Direct:              %.1f ft/s\n" + ...
    "Direct Simmons:      %.1f ft/s\n",V_f_bennett_raw,V_f_simmons_bennett_raw);
fprintf( ...
    "Prandtl-Glauert:     %.1f ft/s\n" + ...
    "P-G Simmons:         %.1f ft/s\n",V_f_bennett_corr,V_f_simmons_bennett_corr);
fprintf("=====%% Difference from TN4197 Direct=====\n");
fprintf( ...
    "TN4197 Simmons:      %+.2f%%\n" + ...
    "TN4197 P-G:          %+.2f%%\n" + ...
    "TN4197 P-G Simmons:  %+.2f%%\n" + ...
    "1/k Direct:          %+.2f%%\n" + ...
    "1/k Direct Simmons:  %+.2f%%\n" + ...
    "1/k P-G:             %+.2f%%\n" + ...
    "1/k P-G Simmons:     %+.2f%%\n",simm_TN4197_diff,TN4197_corr_diff,simm_TN4197_corr_diff,invk_raw_diff,simm_invk_raw_diff,invk_corr_diff,simm_invk_corr_diff);
fprintf( ...
    "U vs. g Direct:      %+.2f%%\n" + ...
    "U vs. g Simmons:     %+.2f%%\n",ug_raw_diff,simm_ug_raw_diff);
fprintf( ...
    "U vs. g P-G:         %+.2f%%\n" + ...
    "U vs. g P-G Simmons: %+.2f%%\n",bennett_corr_diff,simm_bennett_corr_diff);
fprintf( ...
    "Bennett Direct:      %+.2f%%\n" + ...
    "Bennett Simmons:     %+.2f%%\n",bennett_raw_diff,simm_bennett_raw_diff);
fprintf( ...
    "Bennett P-G:         %+.2f%%\n" + ...
    "Bennett P-G Simmons: %+.2f%%\n",ug_corr_diff,simm_ug_corr_diff);
fprintf("=====%% Difference from TN4197 Simmons=====\n");
fprintf( ...
    "TN4197 P-G:          %+.2f%%\n" + ...
    "TN4197 P-G Simmons:  %+.2f%%\n" + ...
    "1/k Direct:          %+.2f%%\n" + ...
    "1/k Direct Simmons:  %+.2f%%\n" + ...
    "1/k P-G:             %+.2f%%\n" + ...
    "1/k P-G Simmons:     %+.2f%%\n",diff2_simm_TN4197_corr,diff2_simm_TN4197_simm_corr,diff2_simm_invk_raw,diff2_simm_invk_simm_raw,diff2_simm_invk_corr,diff2_simm_invk_simm_corr);
fprintf( ...
    "U vs. g Direct:      %+.2f%%\n" + ...
    "U vs. g Simmons:     %+.2f%%\n",diff2_simm_ug_raw,diff2_simm_ug_simm_raw);
fprintf( ...
    "U vs. g P-G:         %+.2f%%\n" + ...
    "U vs. g P-G Simmons: %+.2f%%\n",diff2_simm_ug_corr,diff2_simm_ug_simm_corr);
fprintf( ...
    "Bennett Direct:      %+.2f%%\n" + ...
    "Bennett Simmons:     %+.2f%%\n",diff2_simm_bennett_raw,diff2_simm_bennett_simm_raw);
fprintf( ...
    "Bennett P-G:         %+.2f%%\n" + ...
    "Bennett P-G Simmons: %+.2f%%\n",diff2_simm_bennett_corr,diff2_simm_bennett_simm_corr);

%% Function Definitions

function [V_f] = TN4197(G_E, h, t, C_r, C_t, Theta_LE, p, p_0, a)
    %{
    Calculates flutter velocity for a trapezoidal fin based on the formula
    in NACA TN4197. https://ntrs.nasa.gov/citations/19930085030
    Ultimately: V_f = a * sqrt( G_E / ( ((39.3 * AR^3) / ((t/c)^3 * (AR + 2))) * ((lambda + 1) / 2) * (p/p_0) ) )

    Inputs:
    C_r - Root Chord, Inches
    C_t - Tip Chord, Inches
    Theta_LE - Leading-Edge Sweep angle, Degrees
    h - Fin Height, Inches
    t - Fin Thickness, Inches
    a - Speed of Sound in Fluid, Feet / Second
    G_E - Shear Modulus, Pounds / Inches^2
    p - Local Air Pressure, Pounds / Inches^2
    p_0 - Sea Level Air Pressure, Pounds / Inches^2

    Outputs:
    V_f - Flutter Velocity, Feet / Second
    %}
    
    finArea = h * 0.5 * (C_r + C_t); % area of fin, inches^2
    AR = h^2 / finArea; % Fin aspect ratio, unitless
    lambda = C_t / C_r; % Fin taper ratio, unitless
    
    V_f = a * sqrt(G_E / ( ((39.3 * AR^3) / ((t/C_r)^3 * (AR + 2))) * ((lambda + 1) / 2) * (p/p_0)));
end

function [Uf] = TR496TR685(freq_alpha, freq_h, a_h, x_alpha, r_alpha, b, mu, invkstepsize, invkMax, g_h, g_alpha)
    %{
    Calculates flutter velocity based on the sqrt(X) vs 1/k method.
    Originally found in NACA TR496: https://ntrs.nasa.gov/citations/19930090935
    again in NACA TR685: https://ntrs.nasa.gov/citations/19930091762
    and also in Y.C. Fung's "An Introduction to the Theory of
    Aeroelasticity"
    Flutter condition is when the real and imaginary portions of sqrt(X)
    plotted against 1/k cross

    freq_alpha - Natural pitching (pure rotation about E.A.) frequency of fin, rad/s
    freq_h - Natural plunge (pure bending about E.A.) frequency of fin, rad/s
    a_h - location of elastic axis (E.A.) of fin behind mid-chord divided by semichord length, unitless
    x_alpha - location of c.g. behind E.A. as ratio of semichord, unitless
    r_alpha - reduced radius of gyration around E.A. divided by semichord, unitless
    b - semichord length (half of length of fin), feet (can be other length unit, defines velocity as [unit(b)]/s)
    mu - nondimensional mass ratio, unitless
    g_h - plunge damping coefficient, unitless, typ 0.005 (per FinSim)
    g_alpha - pitch damping coefficient, unitless, typ 0.005 (per FinSim)
    invkstepsize - size of difference between discrete points of 1/k, unitless recommended between 0.0001 - 0.000001
    invkMax - maximum value of 1/k to calculate to, unitless, typ 6-14, adjust higher if the parabola is not closed
    %}
    i = sqrt(-1);
    invkrange = [invkstepsize,invkMax];
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
    
    delta_I_A = (g_h + g_alpha) * mu^2 * r_alpha^2 * freq_h.^2 / freq_alpha^2;
    delta_I_B = (mu .* freq_h.^2 ./ freq_alpha.^2 .* ((g_h .* E_R) + E_I)) + (mu .* r_alpha.^2 .* (A_I + (g_alpha .* A_R)));
    delta_I_C = (A_I .* E_R) - (B_R .* D_I) + (A_R .* E_I) - (B_I .* D_R);
    
    X_R1 = (-delta_R_B - sqrt(delta_R_B.^2 - (4 .* delta_R_A .* delta_R_C))) ./ (2 .* delta_R_A);
    X_R2 = (-delta_R_B + sqrt(delta_R_B.^2 - (4 .* delta_R_A .* delta_R_C))) ./ (2 .* delta_R_A);
    
    if(~delta_I_A == 0)
        X_I1 = (-delta_I_B - sqrt(delta_I_B.^2 - (4 .* delta_I_A .* delta_I_C))) ./ (2 .* delta_I_A);
        X_I2 = (-delta_I_B + sqrt(delta_I_B.^2 - (4 .* delta_I_A .* delta_I_C))) ./ (2 .* delta_I_A);
    else
        X_I1 = -delta_I_C ./ delta_I_B;
        X_I2 = -delta_I_C ./ delta_I_B;
    end
    
    iscomplex = (imag(X_R1) ~= 0);
    X_R1(iscomplex) = NaN;
    iscomplex = (imag(X_R2) ~= 0);
    X_R2(iscomplex) = NaN;
    
    rt_X_R1 = sqrt(X_R1);
    rt_X_R2 = sqrt(X_R2);
    rt_X_I1 = sqrt(X_I1);
    rt_X_I2 = sqrt(X_I2);
    
    solutionmatrix = [k; k_inv; X_R1; X_R2; X_I1; rt_X_R1; rt_X_R2; rt_X_I1; rt_X_I2];
    
    XRatio1 = abs(1 - abs(solutionmatrix(6,:) ./ solutionmatrix(8,:)));
    XRatio2 = abs(1 - abs(solutionmatrix(7,:) ./ solutionmatrix(8,:)));
    
    [~,idx1] = find(XRatio1 == min(XRatio1,[],"omitnan"));
    [~,idx2] = find(XRatio2 == min(XRatio2,[],"omitnan"));
    r1 = XRatio1(:,idx1);
    r2 = XRatio2(:,idx2);
    if(r1 < r2)
        flutterPoint = solutionmatrix(:,idx1);
        flutter_rtX = flutterPoint(6);
        Uf = freq_alpha * b * flutterPoint(2) / flutterPoint(6);
        flutterFreq = freq_alpha / flutterPoint(6) / (2 * pi());
    else
        flutterPoint = solutionmatrix(:,idx2);
        flutter_rtX = flutterPoint(7);
        Uf = freq_alpha * b * flutterPoint(2) / flutterPoint(7);
        flutterFreq = freq_alpha / flutterPoint(7) / (2 * pi());
    end
    inv_kf = flutterPoint(2);
    kf = flutterPoint(1);
end

function [Uf, Ud] = UvsgNoDamp(freq_alpha, freq_h, a_h, x_alpha, r_alpha, b, mu, kStepSize,kMax)
    %{
    freq_alpha - Natural pitching (pure rotation about E.A.) frequency of fin, rad/s
    freq_h - Natural plunge (pure bending about E.A.) frequency of fin, rad/s
    a_h - location of elastic axis (E.A.) of fin behind mid-chord divided by semichord length, unitless
    x_alpha - location of c.g. behind E.A. as ratio of semichord, unitless
    r_alpha - reduced radius of gyration around E.A. divided by semichord, unitless
    b - semichord length (half of length of fin), feet (can be other length unit, defines velocity as [unit(b)]/s)
    mu - nondimensional mass ratio, unitless
    kstepsize - size of difference between discrete points of k, unitless recommended between 0.0001 - 0.000001
    kMax - maximum value of k to calculate to, unitless, typ 6-14, adjust higher if the parabola is not closed
    %}
    i = sqrt(-1);
    kRange = [kStepSize,kMax];
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
end

function [V_f] = SimmonsFlutter(x_bar, r_bar, freq_h, freq_alpha, b, mu, Mach)
   
    % Only valid for supersonic flight based on supersonic lift assumptions
    % r_bar = radius of gyration wrt axis of rotation/elastic axis
    % x_bar = distance of cg from elastic/axis/axis of rotation
    V_f = b * freq_alpha * sqrt((mu * r_bar^2 * sqrt(Mach^2 - 1) / (x_bar * b)) * (1 - (freq_h / freq_alpha)^2)^2 + (4 * (x_bar / r_bar)^2 * (freq_h / freq_alpha)^2));

end

function [V_f] = BennetFlutter(C_r,C_t,Theta_LE,h,t,a,G_E,p,p_0,gamma)
    %{
    Calculates Bending-Torsion Flutter Velocity for a trapezoidal fin based on Bennet's formula (and derivations)
    uploaded to Richard Nakka's rocketry website (Dec 2023)
    Link: https://www.nakka-rocketry.net/articles/Calculating_Fin_Flutter_Velocity_Bennett-12-23.pdf
    Ultimately: V_f = a * sqrt(G_E / (((DN * AR^3) / ((t/C_r)^3 * (AR + 2))) * ((lambda + 1)/2) * (p/p_0)))
    
    Inputs:
    C_r - Root Chord, Inches
    C_t - Tip Chord, Inches
    Theta_LE - Leading-Edge Sweep angle, Degrees
    h - Fin Height, Inches
    t - Fin Thickness, Inches
    a - Speed of Sound in Fluid, Feet / Second
    G_E - Shear Modulus, Pounds / Inches^2
    p - Local Air Pressure, Pounds / Inches^2
    p_0 - Sea Level Air Pressure, Pounds / Inches^2
    Gamma - Ratio of Speciifc Heats of Fluid, Unitless

    Outputs:
    V_f - Flutter Velocity, Feet / Second

    Shortcomings:
    1) Does not account for mass of the fin
    2) Does not account for air mass
    3) *Assumes constant fin thickness
    4) Only validated to Mach 1.5
    5) Assumes motion is Bending-Torsion
    6) Only directly accounts for shear stiffness
    7) Does not account for position of fin C.G. along fin height
    %}
    
    finArea = h * 0.5 * (C_r + C_t); % area of fin, inches^2
    AR = h^2 / finArea; % Fin aspect ratio, unitless
    lambda = C_t / C_r; % Fin taper ratio, unitless
    sweepLength = tand(Theta_LE) * h; % length from root LE to tip LE
    C_x = ((2 * C_t * sweepLength) + (C_t^2) + (sweepLength * C_r) + (C_t * C_r) + (C_r^2)) / (3 * (C_t + C_r)); % Centroid of fin, effectively center of mass location from elading edge, inches
    epsilon = (C_x / C_r) - 0.25; % distance of fin center of mass behind fin quarter-chord line, unitless
    DN = (24 * epsilon * gamma * p_0) / pi(); % "Denominator Constant", Pounds / Inches^2
    V_f = a * sqrt(G_E / (((DN * AR^3) / ((t/C_r)^3 * (AR + 2))) * ((lambda + 1) / 2) * (p / p_0))); % Flutter Velocity, Feet / Second
end