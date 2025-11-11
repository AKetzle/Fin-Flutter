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

%% 3 - Kholodar Nondimensional Transonic Flutter
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

%% 5 - NACA TR496 and TR685 - REDO, USE WEISSHAAR TEXTBOOK MATRIX VERSION

function [V_f] = TR496TR685(freq_alpha, freq_h, r_alpha, r_h, b, c, rho, m, afin)
    %{
    Calculates flutter velocity for a trapezoidal fin based on the formulas
    in NACA TR496 and NACA TR685 for case 1.
    TR496 Link: https://ntrs.nasa.gov/citations/19930090935
    TR685 Link: https://ntrs.nasa.gov/citations/19930091762
    
    inputs:
    rho = density of air in slugs / ft^3
    m = mass of fin in slugs
    r_alpha = reduced radius of gyration for pitch
    afin = location of axis of rotation (elastic axis) of fin from leading
    edge (ratio/coordinate) a = -2(stiff axis)/100 - 1
    x_alpha = distance between fin CG and afin location
    c = distance axis of rotation of control surface from midchord (unused)
    * beta inputs are omitted because this assumes a fin without ailerons
    %}
    r_beta = 0;
    x_beta = 0;
    g_alpha = 0;
    g_beta = 0;
    g_h = 0;
    kappa = pi() * b^2 * rho / m; % mass ratio - mu is typical nowadays but this is to be consistent with the paper so I don't go nuts.

    % T equations from TR496
    p = -(1/3) * (sqrt(1 - c^2))^3;
    T_1 = -(1/3) * sqrt(1 - c^2) * (2 + c^2) + (c * acos(c));
    T_2 = c * (1 - c^2) - sqrt(1 - c^2) * (1 + c^2) * acos(c) + c*acos(c)^2;
    T_3 = -((1/8) + c^2) * acos(c)^2 + (1/4) * c * sqrt(1 - c^2) * acos(c) * (7 + (2 * c^2)) - ((1/8) * (1 - c^2) * (5 * c^2 + 4));
    T_4 = -acos(c) + c * sqrt(1- c^2);
    T_5 = -(1 - c^2) - acos(c)^2 + (2 * c * sqrt(1 - c^2) * acos(c));
    T_6 = T_2;
    T_7 = -((1/8) + c^2) * acos(c) + ((1/8) * c * sqrt(1 - c^2) * (7 + 2 * c^2));
    T_8 = -((1/3) * sqrt(1 - c^2) * (2 * c^2 + 1)) + c * acos(c);
    T_9 = (1/2) * (((1/3) * sqrt(1 - c^2)) + (afin * T_4)); % WARNING: there's an inconsistency in TR496 about this one, look into it!!!
    T_10 = sqrt(1 - c^2) + acos(c);
    T_11 = acos(c)* (1 - 2 * c) + (sqrt(1 - c^2) * (2 - c));
    T_12 = sqrt(1 - c^2) * (2 + c) - (acos(c) * (2 * c + 1));
    T_13 = (1/2) * (-T_7 - ((c - afin) * T_1));
    T_14 = (1/16) + (afin * c / 2);

    A_alpha1 = (r_alpha^2 / kappa) + ((1/8) + afin^2);
    A_alpha2 = 0.5 - afin;
    A_beta1 = (r_beta^2 / kappa) - ((T_7 / pi())) + ((c - afin) * ((x_beta / kappa) - (T_1 / pi())));
    A_beta2 = (1 / pi()) * (-2 * p - ((1/2) - afin) * T_1);
    A_beta3 = (T_4 + T_10) / pi();
    A_h1 = (x_alpha / kappa) - afin;
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

    % Solving Matrices
    
    A_1 = det([A_alpha1, A_h1; C_alpha1, C_h1]);
    B_1 = det([A_alpha1, -(0.5 + afin); C_alpha1, 1]) + ((0.5 - afin) * det([-(0.5 + afin), A_h1; 1, C_h1]));
    C_1 = det([(A_h1 - A_alpha2), -(0.5 + afin); (C_h1 - C_alpha2), 1]);
    D_1 = -det([A_alpha2, A_h1; C_alpha2, C_h1]);

    syms k X;
    for i=1:iterations
        fprintf("Iteration: %d\n     X0: %f\n     k0: %f\n",i,X0,k0);
        % Bessel Function Stuff - TR496 cites Niels Nielsen's 1904 book about
        % these functions, using the matlab built-in with the example in TR685
        % the calculated values line up.
        J_1 = besselj(1,k);
        J_0 = besselj(0,k);
        Y_0 = bessely(0,k);
        Y_1 = bessely(1,k);
        
        F = (J_1 * (J_1 + Y_0) + (Y_1 * (Y_1 - J_0))) / ((J_1 +Y_0)^2 + (Y_1 - J_0)^2);
        G = - ((Y_1 * Y_0) + (J_1 * J_0)) / ((J_1 + Y_0)^2 + (Y_1 - J_0)^2);
        
        % I and R functions
        
        R_aalpha = -A_alpha1 + ((0.25 - afin^2) * (2 * G / k)) - ((0.5 + afin) * (2 * F / k^2));
        %R_abeta = -A_beta1 + (A_beta3 / k^2) + (0.5 + a) * ((T_11 * G / (pi() * k)) - (T_10 * F * 2 / (pi() * k^2)));
        %R_ah = -A_h1 + (0.5 + a) * (2 * G / k);
        
        %R_balpha = -B_alpha1 - (T_12 * (((0.5 - a) * 2 * G / k) - (2 * F / k^2)) / (2 * pi()));
        %R_bbeta = -B_beta1 + (B_beta3 / k^2) - (T_12 * ((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / (pi() * k^2))) / (2 * pi()));
        %R_bh = -B_h1 - (T_12 * G / k) / pi();
        
        %R_calpha = -C_alpha1 - ((0.5 - a) * 2 * G / k) + (2 * F / k^2);
        %R_cbeta = -C_beta1 - (T_11 * G / (pi() * k)) + (T_10 * 2 * F / (pi() * k^2));
        R_ch = -C_h1 - (2 * G / k);
        
        I_aalpha = (1/k) * (-(((1/2) + afin) * 2 * G / k) - (((1/4) - afin^2) * 2 * F) + A_alpha2);
        %I_abeta = (-(0.5 + a) * (((T_10 * G * 2) / (pi() * k)) + (T_11 * F / pi())) + A_beta2) / k;
        %I_ah = -(0.5 + a) * 2 * F;
        
        %I_balpha = ((T_12 / (2 * pi())) * ((2 * G / k) + ((0.5 - a) * (2 * F))) + B_alpha2) / k;
        %I_bbeta = (((T_12 / (2 * pi())) * ((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / pi()))) + B_beta2) / k;
        %I_bh = (T_12 * F / pi()) / k;
        
        %I_calpha = ((2 * G / k) - ((0.5 - a) * 2 * F) + C_alpha2) / k;
        %I_cbeta = (((T_11 * 2 * G / (2 * pi() * k)) - (T_10 * 2 * F / pi())) + C_beta2) / k;
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
    
        if(Xdiff < convergence & kdiff < convergence)
            fprintf("Converged in %d iterations\n",i);
            break
        elseif(abs(Xdiff) > 10 & abs(kdiff) > 10)
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
    V_f = r_alpha * freq_alpha * b / (sqrt(kappa) * k * sqrt(X));
    Mi = V_f / amach;
    V_f_c = amach * Mi * sqrt(sqrt(Mi^4 + 4) - Mi^2) / sqrt(2); % compressibility correction
end

%% 6 - Simmons simplification of Kearns (@ John Hopkins) equation

function [V_f] = SimmonsFlutter(x_bar, r_bar, freq_h, freq_alpha)
   
    % Only valid for supersonic flight based on supersonic lift assumptions
    % r_bar = radius of gyration wrt axis of rotation/elastic axis
    % x_bar = distance of cg from elastic/axis/axis of rotation
    V_f = b * freq_alpha * sqrt((mu * r_bar^2 * sqrt(Mach^2 - 1) / (x_bar * b)) * (1 - (freq_h / freq_alpha)^2)^2 + (4 * (x_bar / r_bar)^2 * (freq_h / freq_alpha)^2));

end