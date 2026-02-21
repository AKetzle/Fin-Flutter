clc, clear;

syms Mc Mi

eq1 = Mc^2 / sqrt(Mc^2 - 1) == Mi^2;
assume(Mc >= 1)
assumeAlso(Mc,"real")
assumeAlso(Mi,"real")
assumeAlso(Mi,"positive")
assumeAlso(Mi >= 0)
correction_eq1 = solve(eq1,Mc,"Real",true)

eq2 = Mc^2 / sqrt(1 - Mc^2) == Mi^2;
assume(Mc <= 1)
assumeAlso(Mc >= 0)
assumeAlso(Mc,"real")
assumeAlso(Mi,"real")
assumeAlso(Mi,"positive")
correction_eq2 = solve(eq2,Mc,"Real",true)
