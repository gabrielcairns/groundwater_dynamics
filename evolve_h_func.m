function fun = evolve_h_func(hUV, h_current, Ps, b, S, xg, xgt, dX, J, delta_rho, K, dt)
% Finds h (and complementarity variables U and V) at the next timestep
% Ps, b, S, xg, xgt are ice overburden, basement, upper surface, GL, rate
% of change of GL at the future (sought) timestep
% dX is grid spacing, J is number of mesh points (this currently assumes a uniform grid)
% delta and Kay are parameters
% dt is the time step

% Complementarity variables
h = hUV(1:J+1);
U = hUV(J+2:2*J+1);   V = hUV(2*J+2:3*J+1); 

% Reconstruction of mesh grid
XH = dX*((1:J)-1/2)'; % X at the the halfway grid points

% Midpoint h
h_mid = (h(1:J)+h(2:J+1))/2;

% Hydraulic gradient
PiX = (Ps(2:J+1)-Ps(1:J)+S(2:J+1)-S(1:J))/dX;

% Advection speed
adv = K*(PiX/(xg^2) + delta_rho*(b(2:J+1)-b(1:J))/(dX*xg^2) ) + (xgt/xg)*XH ;

% Upwind value of h
h_upw = h(1:J).*(adv<=0)+h(2:J+1).*(adv>0);

% Flux into adjacent cells (midpoint diffusion for second term)
salt_flux = [0; h_upw.*adv + K*delta_rho*h_mid.*(h(2:J+1)-h(1:J))/(dX*xg^2)] ;

% Right hand side of equation
DhDT =  (salt_flux(2:J+1)-salt_flux(1:J))/dX - (xgt/xg)*h(1:J)  ;

% Calculated h, but not at the final meshpoint
fun(1:J) = U.*V; % Solving the actual problem
fun(J+1) = h(J+1) + b(J+1) - S(J+1); % BC at X=1
fun(J+2:2*J+1) =  S(1:J) - b(1:J) - h(1:J) - U.^2; % Complimentarity 1
fun(2*J+2:3*J+1) =  h(1:J) - h_current(1:J) - dt*DhDT + V.^2; % Complimentarity 2
end
