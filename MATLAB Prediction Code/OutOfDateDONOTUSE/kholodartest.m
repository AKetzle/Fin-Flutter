clc, clear;

freq_h = 0.8;
freq_alpha = 1;
x_alpha = 0.25;
r_alpha = sqrt(0.75);
a = -0.6;
m = 1;
rho = 0.002378;
b = 6;

KholodarTransonic(freq_h,freq_alpha,r_alpha,x_alpha,0,0,0,0,b,m,rho)

function [V_f] = KholodarTransonic(freq_h, freq_alpha, r_alpha, x_alpha, c_lh_bar, c_lalpha_bar, c_mh_bar, c_malpha_bar, b, m, rho)
    syms freq_bar h_bar alpha_bar V

    mu = pi() * b^2 * rho / m;
    m1 = -freq_bar^2 * [1, x_alpha; x_alpha, r_alpha^2];
    m2 = 4 * [(freq_alpha^2 / freq_h^2), 0; 0, r_alpha^2] / V^2;
    m3 = 4 * [c_lh_bar, c_lalpha_bar; (-2 * c_mh_bar), (-2 * c_malpha_bar)] / (pi() * mu);
    F = (m1 + m2 + m3);
    assume(V > 0);
    assume(freq_bar > 0);
    freq_bar = solve(det(F) == 0 + 0i,freq_bar)
    fluttereq = F * [h_bar; alpha_bar] == [0; 0]

end