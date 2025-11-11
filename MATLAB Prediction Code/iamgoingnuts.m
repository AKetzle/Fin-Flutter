clc, clear;
sympref("FloatingPointOutput",true);
syms X_r X_i k_1 X;
realeq = (0.25 * X_r^2) - (5.15213 * X_r) + 5.76524 == 0;
imageq = (1.28454 * X_i) - 3.30557 == 0;
sqrtxr = sqrt(solve(realeq, X_r))
sqrtxi = sqrt(solve(imageq,X_i))
test = (1.28454 * X) - 3.30557 == (0.25 * X^2) - (5.15213 * X) + 5.76524;
sqrt(solve(test,X)) 