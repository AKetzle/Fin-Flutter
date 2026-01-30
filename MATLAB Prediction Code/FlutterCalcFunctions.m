%% Various available flutter calculation functions
% Alexander Ketzle

%% 1 - Bennett 2023 Flutter Equation (From Richard Nakka's Website)

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

%% 2 - Weisshaar Quasi-Steady Flutter Equation

%% 3 - (NOT FOR USE) Kholodar Nondimensional Transonic Flutter
function [V_f] = KholodarTransonic(freq_h, freq_alpha, r_alpha, c_lh_bar, c_lalpha_bar, c_mh_bar, c_malpha_bar, b, m, rho)
    syms freq_bar h_bar alpha_bar V

    mu = pi() * b^2 * rho / m;
    m1 = -freq_bar^2 * [1, x_alpha; x_alpha, r_alpha^2];
    m2 = 4 * [(freq_alpha^2 / freq_h^2), 0; 0, r_alpha^2] / V^2;
    m3 = 4 * [c_lh_bar, c_lalpha_bar; (-2 * c_mh_bar), (-2 * c_malpha_bar)] / (pi() * mu);
    fluttereq = (m1 + m2 + m3) * [h_bar; alpha_bar] == [0; 0];
    V_f = vpasolve(fluttereq,V);

end
%% 4 - NACA TN4197

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
    
    V_f = a * sqrt( G_E / ( ((39.3 * AR^3) / ((t/C_r)^3 * (AR + 2))) * ((lambda + 1) / 2) * (p/p_0)));
end

%% 5 - Inverse K method - Based on NACA TR496, NACA TR685, and Y.C. Fung's Theory of Aeroelasticity

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

%% 6 - Simmons simplification of Kearns (@ John Hopkins) equation

function [V_f] = SimmonsFlutter(x_bar, r_bar, freq_h, freq_alpha, b, mu, Mach)
   
    % Only valid for supersonic flight based on supersonic lift assumptions
    % r_bar = radius of gyration wrt axis of rotation/elastic axis
    % x_bar = distance of cg from elastic/axis/axis of rotation
    V_f = b * freq_alpha * sqrt((mu * r_bar^2 * sqrt(Mach^2 - 1) / (x_bar * b)) * (1 - (freq_h / freq_alpha)^2)^2 + (4 * (x_bar / r_bar)^2 * (freq_h / freq_alpha)^2));

end

%% 7 - U vs. g method no structural damping - Based on Weisshaar formulation of Theodorsen equations

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