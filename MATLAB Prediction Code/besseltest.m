clc, clear;
syms k
J_1 = besselj(1,1)
J_0 = besselj(0,1)
Y_0 = bessely(0,1)
Y_1 = bessely(1,1)

F = (J_1 * (J_1 + Y_0) + (Y_1 * (Y_1 - J_0))) / ((J_1 +Y_0)^2 + (Y_1 - J_0)^2)
G = - ((Y_1 * Y_0) + (J_1 * J_0)) / ((J_1 + Y_0)^2 + (Y_1 - J_0)^2)