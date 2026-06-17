o2 = 0:5:400;
T = 0:1:30;
S = 35;
[T,S] = ndgrid(T,S);
o2_sol = oxy_sol(T,S);
load('float_corr_Mar-2026_D_A.mat');
o2_4 = 

plot()