% Plot a case of 'pocket' steady states
% Fig. 4

k0 = 1e-12; 
[delta_rho, z_sea, Ld_km, Hd, Td_kyr, r_i, k0_to_Kay, qEd_mmyr] = default_params();
alph = 0.05;
Kay = k0_to_Kay*k0; 

% Meshgrid
J = 1000;      xg = 0.5; 
dx = xg/J;     x = 0:dx:xg ; 
xtra = [x 0.6];

% Basement geometry
b =  -2.5 + exp( -(x-.25).^2/.025^2) ;
bf = @(x) -2.5 + exp( -(x-.25).^2/.025^2 );
H = -1 - b;     S = H + b;

% Ice sheet
ps = r_i*SSA_ice(x, S, alph); 

% Least salty steady state
hs = (ps(end) - ps - H + (1+delta_rho) * (H(end) + b(end) - b) )/delta_rho;
hs_uncut = hs; 
xn = max(x(hs<0));      hs(x<=xn)=0;

figure(1); clf; set(gcf(),"Position",[1 1 775 700]);
ah = init_axes(subplot(4, 1, [1 2]),'','$z$ / m',[0 Ld_km*max(xtra)],Hd*[min(b) max(H+b+ps)+0.1]);
aqE = init_axes(subplot(4, 1, 3),'','$q_E$ / mm yr$^{-1}$',[0 Ld_km*max(xtra)],[-30 30]);   
asl = init_axes(subplot(4, 1, 4),'$x$ / km','$z$ / m', [0 Ld_km*max(xtra)],[-2500 -1000]);
grid(aqE, 'on');   

% Plot lines pS + S + delta z = cst.
c = linspace(-1.5, 3.5, 21);
for C=c
    hc = hs_uncut + b + C;
    hc(hs_uncut + C > H) = NaN;
    p3 = plot(ah, Ld_km*x, Hd*hc, 'Color', [1 0.5 0.5], 'linestyle', '--', 'LineWidth', 0.6);
end

% Basement, ice, ocean
plot(ah, Ld_km*[flip(xtra) xtra], Hd*[flip(b) b(1) H(1)+b(1) H+b], 'k-' , 'linewidth', 2);  
fill(ah, Ld_km*[flip(xtra) xtra], Hd*[flip(b) b(1) min(b).*ones(size(xtra))], [.5 .5 .5]);  
pps = plot(ah, Ld_km*[x 0.575 0.575 x(end)], ...
    Hd*[H+b+ps/r_i (H(end)+b(end))+[0.95 .2 0]*ps(end)/r_i], 'c-', 'linewidth', 2);
plot(ah, Ld_km*[0.575 0.6], [0 0], 'b-', 'linewidth', 2)

% Where can pockets form?
takeOff = ps(2:end) + H(2:end) + (1+delta_rho)*b(2:end)>...
    ps(1:end-1) + H(1:end-1) + (1+delta_rho)*b(1:end-1);
XP = x(takeOff);

for xp = [XP(end:-5:1) XP(1)]
    % Alternative pocket state
    h_alt = hs;
    h_alt(x<xp) = (ps(x==xp)-ps(x<xp) + H(x==xp)-H(x<xp)...
        + (1+delta_rho)*(b(x==xp)-b(x<xp)))/delta_rho;
    h_alt(h_alt<0)=0;

    % Plot the new pocket SS
    if xp==XP(end)
        h_maxp = h_alt;
        pmaxp = plot(ah, Ld_km*x, Hd*(h_alt+b), 'r-', 'LineWidth', 3);
    else
        plot(ah, Ld_km*x, Hd*(h_alt+b), 'r-', 'LineWidth', 1);
    end
end
pbottom = plot(ah, Ld_km*x, Hd*(hs+b), 'r-', 'LineWidth', 3);

% Legend
legend(ah, [pmaxp, pps, p3], {'Steady states', 'Ice surface $S+H_i$' , '$p_S+S+\delta z =$ cst.'}, ...
    'fontsize', 14, 'interpreter', 'latex', 'location', 'northwest');

% velocity, flux, qE 
uf = ( ps(2:J+1)-ps(1:J) + S(2:J+1)-S(1:J) )/dx;
us = uf + delta_rho*( b(2:J+1)-b(1:J) + hs(2:J+1)-hs(1:J) )/dx ;   
usp = uf + delta_rho*( b(2:J+1)-b(1:J) + h_maxp(2:J+1)-h_maxp(1:J) )/dx ; 
G = H(1:J).*uf + 0.5*(hs(2:J+1)+hs(1:J)).*(us-uf);
Gp = H(1:J).*uf + 0.5*(h_maxp(2:J+1)+h_maxp(1:J)).*(usp-uf);
qE = Kay*(G(2:end)-G(1:end-1))/dx;
qEp = Kay*(Gp(2:end)-Gp(1:end-1))/dx;

% Plot (fresh) streamlines psi = cst. given by z - s = psi/u 
psi_max = max(-[0 uf].*(H-hs));
psis = 0:psi_max/10:psi_max; 
for psi = psis
    plot(asl, Ld_km*x(2:end), Hd*(b(2:end)+hs(2:end)-psi./uf), 'b-', 'LineWidth', 2.5);
end
for psi = psis
    plot(asl, Ld_km*x(2:end), Hd*(b(2:end)+h_maxp(2:end)-psi./uf), 'g-', 'LineWidth', 1.5);
end
% Basement geometry, limits and labels
plot(asl, Ld_km*[flip(xtra) xtra], Hd*[flip(b) b(1) H(1)+b(1) H+b], 'k-' , 'linewidth', 2);  
fill(asl, Ld_km*[flip(xtra) xtra], Hd*[flip(b) b(1) min(b).*ones(size(xtra))], 'k', 'facealpha', 0.5);  

plot(aqE, Ld_km*[0 max(xtra)], [0 0], 'color', [.5 .5 .5], 'LineWidth', 1);
pm = plot(aqE, Ld_km*x(2:end-1), qEd_mmyr*qE, 'b-', 'LineWidth', 3);
pma = plot(aqE, Ld_km*x(2:end-1), qEd_mmyr*qEp, 'g-', 'LineWidth', 1.5);
legend(asl, [pm pma], {'No pocket', 'Maximal pocket'}, 'fontsize', 14 , 'interpreter', 'latex', 'location', 'southeast');

% Subplot labels
text(ah, 5, -1250, '(a)', 'fontsize', 14, 'backgroundcolor', 'w', 'interpreter', 'latex');
text(aqE, 5, 20, '(b)', 'fontsize', 14, 'interpreter', 'latex');
text(asl, 5, -1300, '(c)', 'fontsize', 14, 'backgroundcolor', 'w', 'interpreter', 'latex');

asl.YLabel.Position(1) = -35;
ah.YLabel.Position(1) = -35;
aqE.YLabel.Position(1) = -35;