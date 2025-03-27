% o2satv2b.m                                       by:  Edward T Peltzer, MBARI
%                                                  revised:  2013 Sep 23.
%
% CALCULATE OXYGEN CONCENTRATION AT SATURATION
%
% Source:  Garcia & Gordon (1992).  Oxygen solubility in seawater:
%          Better fitting equations.  L&O 37: 1307-1312.
%
% Input:       S = Salinity (pss-78)
%              T = Temp (deg C)
%
% Output:      Oxygen saturation at one atmosphere (umol/kg).
%
%                        O2 = o2satv2b(S,T).

function [O2] = o2satv2b(S,T)


% DEFINE CONSTANTS, ETC FOR SATURATION CALCULATION

%  The constants used are for units of umol O2/kg.

  A0 =  5.80818;
  A1 =  3.20684;
  A2 =  4.11890;
  A3 =  4.93845;
  A4 =  1.01567;
  A5 =  1.41575;
  
  B0 = -7.01211e-03;
  B1 = -7.25958e-03;
  B2 = -7.93334e-03;
  B3 = -5.54491e-03;
  
  C0 = -1.32412e-07;
  
%   Calculate Ts from T (deg C)  

  Ts = log((298.15 - T) ./ (273.15 + T));
  
%   Calculate O2 saturation in umol O2/kg.

  A = ((((A5 .* Ts + A4) .* Ts + A3) .* Ts + A2) .* Ts + A1) .* Ts + A0;
  
  B = ((B3 .* Ts + B2) .* Ts + B1) .* Ts + B0;
  
  O2 = exp(A + S.*(B + S.*C0));
