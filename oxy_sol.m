function o2_sol = oxy_sol(t,s)
Ts = log((298.15 - t)./(273.15 + t)); 
% The coefficents below are from the second 
% column of Table 1 of Garcia and Gordon (1992)
a0 =  5.80871; 
a1 =  3.20291;
a2 =  4.17887;
a3 =  5.10006;
a4 = -9.86643e-2;
a5 =  3.80369;
b0 = -7.01577e-3;
b1 = -7.70028e-3;
b2 = -1.13864e-2;
b3 = -9.51519e-3;
c0 = -2.75915e-7;
o2_sol = exp(a0 + Ts.*(a1 + Ts.*(a2 + Ts.*(a3 + Ts.*(a4 + a5*Ts)))) + ...
     s.*(b0 + Ts.*(b1 + Ts.*(b2 + b3*Ts)) + c0*s));