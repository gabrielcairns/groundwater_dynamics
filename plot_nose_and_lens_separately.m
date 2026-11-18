% Creates plots of nose and lens steady states (inc. qE, streamlines)
% Fig. 2

k0 = 1e-12; 
[delta_rho, zsea, Ldkm, Hd, Tkyr, ri, k0_to_Kay, qEdmmyr] = default_params();
Kay = k0_to_Kay*k0; 

% Meshgrid
J = 1000;   xg = 0.4;   dx = xg/J;     
xnose = 0:dx:xg ; % has length J+1
xlens = xnose(xnose<=0.2);
xtral = [xlens 0.25];
xtran = [xnose 0.5]; % for plotting shelf / etc

% Basement geometry
bl = -3*ones(size(xlens));     Hl = 2*ones(size(xlens)); 
bn = -3*ones(size(xnose));     Hn = 2*ones(size(xnose)); 
Sl = Hl+bl;     Sn = Hn+bn;    

% Ice sheets 
psl = ri*SSA_ice(xlens, Hl+bl, .1);
psn = ri*SSA_ice(xnose, Hn+bn, .1);

% Steady states
SSl = max((psl(end) - psl - Hl + (1+delta_rho)*(Hl(end) + bl(end) - bl) )/delta_rho, 0);
SSn = max((psn(end) - psn - Hn +(1+delta_rho)*(Hn(end) + bn(end) - bn) )/delta_rho, 0);

% Velocities, fluxes, qE
ufl = -( psl(2:end)-psl(1:end-1) + Sl(2:end)-Sl(1:end-1) )/dx;
ufn = -( psn(2:J+1)-psn(1:J) + Sn(2:J+1)-Sn(1:J) )/dx;
usl = ufl - delta_rho*(bl(2:end)-bl(1:end-1) + SSl(2:end)-SSl(1:end-1))/dx ; 
usn = ufn - delta_rho*(bn(2:J+1)-bn(1:J) + SSn(2:J+1)-SSn(1:J))/dx ;        
Gl = - Hl(1:end-1).*ufl - 0.5*(SSl(2:end)+SSl(1:end-1)).*(usl - ufl);
Gn = - Hn(1:J).*ufn - 0.5*(SSn(2:J+1)+SSn(1:J)).*(usn - ufn);
qEl = Kay*(Gl(2:end) - Gl(1:end-1))/dx;
qEn = Kay*(Gn(2:end) - Gn(1:end-1))/dx;

% Create figure and axes
figure(1);    clf; set(gcf(),'Position',[1 1 1500 400])
ahl = init_axes(subplot(3, 2, [1 3]),'','$z$ / m',Ldkm*[0 max(xtral)]);
ahn = init_axes(subplot(3, 2, [2 4]),'','$z$ / m',Ldkm*[0 max(xtran)]);
aqel = init_axes(subplot(3, 2, 5),'$x$ / km','$q_E$ / mm yr$^{-1}$',Ldkm*[0 max(xtral)], [-30 15]);
aqen = init_axes(subplot(3, 2, 6),'$x$ / km','$q_E$ / mm yr$^{-1}$', Ldkm*[0 max(xtran)], [-30 15]);

% Basin outline
plot(ahl, Ldkm*[flip(xtral) xtral], Hd*[flip(bl) bl(1) Sl(1) Sl], 'k-', 'linewidth', 2); 
plot(ahn, Ldkm*[flip(xtran) xtran], Hd*[flip(bn) bn(1) Sn(1) Sn], 'k-', 'linewidth', 2); 
% Ice sheets, shelves, ocean
pps = plot(ahl, Ldkm*xlens, Hd*(Hl+bl+psl/ri), 'c-', 'linewidth', 2);
plot(ahn, Ldkm*xnose, Hd*(Hn+bn+psn/ri), 'c-', 'linewidth', 2);
plot(ahl, Ldkm*(xlens(end)+0.5*[0 0.07 0.07 0]), ...
    Hd*(Hl(end)+bl(end)+[0 0.2 0.95 1]*psl(end)/ri), 'c-', 'linewidth', 2);
plot(ahn, Ldkm*(xnose(end)+[0 0.07 0.07 0]), ...
    Hd*(Hn(end)+bn(end)+[0 0.2 0.95 1]*psn(end)/ri), 'c-', 'linewidth', 2);
plot(ahl, Ldkm*[xlens(end)+0.035 xtral(end)], [0 0], 'b-', 'linewidth', 2 );
plot(ahn, Ldkm*[xnose(end)+0.07 xtran(end)], [0 0], 'b-', 'linewidth', 2 );

% Streamlines psi = cst. given by z - s = psi/u
zsll = max([0 ufl].*(Hl-SSl))*(0:0.2:1)'./[0 ufl] ;
zsln = max([0 ufn].*(Hn-SSn))*(0:0.2:1)'./[0 ufn] ;
zsll(SSl+zsll > Hl) = NaN;
zsln(SSn+zsln > Hn) = NaN;

% Plot streamlines
psl = plot(ahl, Ldkm*xlens, Hd*(bl+SSl+zsll), 'Color', .4*[1 1 1], 'LineWidth', 1.5);
plot(ahn, Ldkm*xnose, Hd*(bn+SSn+zsln), 'Color', .4*[1 1 1], 'LineWidth', 1.5);

% Plot steady state s(x, t)
pS = plot(ahl, Ldkm*xlens, Hd*(SSl+bl), 'r-', 'LineWidth', 3);
plot(ahn, Ldkm*xnose, Hd*(SSn+bn), 'r-', 'LineWidth', 3);
% Mark xg1, xg2
scatter(ahl, Ldkm*xlens(end), Hd*(bl(end)+Hl(end)), 32, 'filled', 'k');
scatter(ahn, Ldkm*xnose(end), Hd*(bn(end)+Hn(end)), 32, 'filled', 'k');
text(ahl, Ldkm*xlens(end), Hd*(bl(end)+Hl(end)+0.25), '$x_{g1}$', 'fontsize', 14, 'interpreter', 'latex');
text(ahn, Ldkm*xnose(end), Hd*(bn(end)+Hn(end)+0.25), '$x_{g2}$', 'fontsize', 14, 'interpreter', 'latex');

% Legends and labels
legend(ahl, [pS, pps, psl(1)], {'Steady states $h(x)$', 'Ice surface $S(x)+H_i(x)$', 'Streamlines'}, 'fontsize', 14, ...
    'interpreter', 'latex', 'location', 'southeast');

% Plot qE
plot(aqel, Ldkm*[0 max(xtral)], [0 0], 'color', .4*[1 1 1], 'LineWidth', 1);
plot(aqen, Ldkm*[0 max(xtran)], [0 0], 'color', .4*[1 1 1], 'LineWidth', 1);
plot(aqel, Ldkm*xlens(2:end-1), qEdmmyr*qEl, 'b-', 'LineWidth', 2.5);
plot(aqen, Ldkm*xnose(2:end-1), qEdmmyr*qEn, 'b-', 'LineWidth', 2.5);

% Subplot labels
text(ahl, 5, -200, '(a)(i)', 'fontsize', 14, 'backgroundcolor', 'w', 'interpreter', 'latex');
text(aqel, 5, 5, '(b)(i)', 'fontsize', 14, 'interpreter', 'latex');
text(ahn, 5, -200, '(a)(ii)', 'fontsize', 14, 'backgroundcolor', 'w', 'interpreter', 'latex');
text(aqen, 5, 5, '(b)(ii)', 'fontsize', 14, 'interpreter', 'latex');

% Align the y axis labels with each other
ahl.YLabel.Position(1) = -20;
aqel.YLabel.Position(1) = -20;
ahn.YLabel.Position(1) = -40;
aqen.YLabel.Position(1) = -40;
aqel.XLabel.Position(2) = -40.5;
aqen.XLabel.Position(2) = -40.5;