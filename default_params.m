function [delta_rho, z_sea, Ld_km, Hd, Td_kyr, r_i, k0_to_Kay, qEd_mmyr] = default_params()
% Default parameter values

delta_rho = 1/40; 
phi = 0.3;
Ld = 5e5;   
Hd = 1e3;   
Ld_km=Ld/1e3; 
z_sea=0;
r_i = 0.917;
qEd_mmyr = 10;

Td_kyr = 100;
k0_to_Kay = (Td_kyr * 3.15e10 * 1e3 * 9.8 * Hd)/(1e-3 * phi * Ld^2);
%               [t in s]      rho   g    [z]   k   / mu    phi  [x]^2
end