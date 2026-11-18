% Plot the results of ross_sea_basin.mat
% Fig. 9-11
% Uses cmocean colormaps

[delta_rho, z_sea, Ld_km, Hd, Td_kyr, r_i, k0_to_Kay, qEd_mmyr] = default_params();

load('bHdata.mat', 'St_avg', 'bt_avg','xi_transect');
x_slw = 425/Ld_km;  xmax = max(xi_transect)/Ld_km; 

b = @(x) interp1(xi_transect/Ld_km, bt_avg/Hd, x); 
S = @(x) interp1(xi_transect/Ld_km, St_avg/Hd, x);
H = @(x) S(x) -b(x);

% Results are re-interpolated
JS = 100; % # Interpolation points for s
JE = 100; % # Interpolation points for qE
K0 = 201; % # Interpolation points for timeseries

j_spacetime = 2; % j for which we make space-time plot of qE

% Values of k
%ks = {'1e-10', '3e-12', '1e-12', '3e-13', '1e-13', '3e-14'};
ks = [1e-10 3e-12 1e-12 3e-13 1e-13 3e-14]';
Kays = zeros(6,1);

% qE and s at SLW site
qESLW = zeros(K0,6);
sSLW = zeros(K0,6);

hPLT = zeros(JS,4,4);
xavgqE = zeros(K0,6);

% For each value of k, load and reformat the results
for j= 1:6

    load(sprintf('ross_results/ross_%g.mat', ks(j)),'Kay','hplot','kgif','xg','J','X','h_curr','qE','t','dX','kplots','PS');
    Kays(j) = Kay;
    if j<=4 
        % Extract h data at certain steps (for larger K)
        hPLT(:,:,j) = interp1(hplot(:,[1 2 4 6]),linspace(1,length(hplot),JS)); 
    end

    % Measure at SLW
    sslw = nan(size(kgif)); 
    qEslw = nan(size(kgif)); 
    for k=1:length(kgif)
        xgk = xg(kgif(k));
        if xgk>x_slw
            sslw(k) = b(x_slw)+interp1(xgk*X,h_curr,x_slw);
            qEslw(k) = interp1(xgk*X(2:J),qE(:,k),x_slw);  
        end      

    end
     
    % Creation of mesh for space-time plot
    if j == j_spacetime
        [tmesh, Xmesh] = meshgrid(t(kgif), X(2:end-1));
        xmesh = Xmesh.*xg(kgif);
        qE_st = qE;
    end
    
    % Scale and reinterpolate to time mesh
    xavgqE(:,j) = qEd_mmyr*interp1(t(kgif),sum(qE.*xg(kgif)*dX,1)./xg(kgif),linspace(-1,0,K0));
    qESLW(:,j) = qEd_mmyr*interp1(t(kgif),qEslw,linspace(-1,0,K0));
    sSLW(:,j) = Hd*interp1(t(kgif),sslw,linspace(-1,0,K0));
end

% Depth of lens (m)
disp(Hd*S(x_slw)-sSLW(end,:));

% ross_sea_basin plots at 6 points, we only use 4
kplots = kplots([1 2 4 6]);

% Create new X mesh
X = linspace(0, max(X), JS); 
XE = linspace(0, max(X), JE);  
Kays = round(Kays, 3, 'significant');

% Colours used to plot different values of k
clr = cmocean('matter',7);
clr = clr(2:7,:);
ls = {'-.','-','-.','-','-.','-'};
letters = 'abcd';

% Legend text
ktext = cell(6,1);
for i = 1:6
    ktext{i} = sprintf('$k=%g$ m$^2$, $K=%g$',ks(i), Kays(i));
end

% Snapshots of s(x,t) at points in time
figure(1); clf; set(gcf(),"Position",[1 1 750 750])
for P = 1:length(kplots)
    xgkp = xg(kplots(P));

    ah = init_axes(subplot('Position',[0.12 0.93-.19*P 0.85 0.15]),'','$z$ / m',[0 Ld_km*xmax],[-5400 2600]);
    title(ah,sprintf('(%s) $t=%g$ kyr before present',letters(P),-round(Td_kyr*t(kplots(P)))),'FontSize',16,'Interpreter','latex')

    % SLW position
    plot(ah, Ld_km*x_slw*[1 1], [-5400 2300], 'k--','linewidth',1);
    % Plot outline of sedimentary basin
    plot(ah, Ld_km*xmax*X, Hd*b(xmax*X), 'k-' ,'linewidth', 1.5);
    plot(ah, Ld_km*xmax*X, Hd*S(xmax*X), 'k-','linewidth', 1.5);
   
    % Ice
    plot(ah, Ld_km*xgkp*X, Hd*(S(xgkp*X)+ PS(X,xgkp)/r_i), 'c-','linewidth',2.5);
    shl = 0.03; % Shelf length
    plot(ah, Ld_km*(xgkp+ shl*[0 1 1 0]), ...
        Hd*(H(xgkp)+b(xgkp)+[0 0.2 0.95 1]*(1+delta_rho)*(-H(xgkp)-b(xgkp))/r_i), 'c-','linewidth',2.5);
    plot(ah, Ld_km*[xgkp+shl xmax], [0 0], 'b-','linewidth',2.5);
    
    % Interface(s)
    for j=1:4
        plot(ah, Ld_km*xgkp*X, Hd*(hPLT(:,P,j)'+b(xgkp*X)),'Color',clr(j,:),'LineStyle',ls{j},'linewidth', 1);
    end

    if P==4
        xlabel(ah, '$x$ / km','Fontsize',14,'interpreter','latex');
    else
        ah.XTickLabels = '';
    end

    ah.YLabel.Position(1) = -150;
end

% Legend
aleg = subplot('Position',[0.12 0.04 0.85 0.1]);    aleg.Visible = 'off';
legend(aleg,[ah.Children(4) ah.Children(3) ah.Children(2) ah.Children(1)], ...
    ktext{1:4},'interpreter','latex','Position',[0.18 0.01 0.7 0.08],'numcolumns',2,'fontsize',16);

% Timeseries of average q_e
figure(4);  clf; set(gcf(),'Position',[1 1 750 650])
axg = init_axes(subplot(4,1,1),'','$x_g$ / km'); 
aqeK = init_axes(subplot(4,1,2),'','$\bar{q}_E / K$','auto',[-8 8]); 
aqe = init_axes(subplot(4,1,[3 4]),'$t$ / kyr before present','$\bar{q}_E$ / mm yr$^{-1}$','auto',[-60 60]); 
grid(axg,'on');   
grid(aqe,'minor');  grid(aqeK,'minor');

% Grounding line plot
plot(axg, [-100 Td_kyr*t(kplots)],Ld_km*[xg(1) xg(kplots)], 'k-','linewidth',2);
scatter(axg, Td_kyr*t(kplots),Ld_km*xg(kplots), 32,'k','filled');

% Plot average q_e
for j = 2:6
    plot(aqe, Td_kyr*linspace(-1,0,K0),xavgqE(:,j), 'Color',clr(j,:),'LineStyle',ls{j},'linewidth',1); 
end
for j = 1:6
    plot(aqeK, Td_kyr*linspace(-1,0,K0),.1*xavgqE(:,j)./Kays(j), 'Color',clr(j,:),'LineStyle',ls{j},'linewidth',1);
end

% Labels and legends
legend(aqe,aqeK.Children(6:-1:1), ktext,'location','southoutside','FontSize',14,'Interpreter','latex','NumColumns',2);

aqe.YLabel.Position(1) = -110;
axg.YLabel.Position(1) = -110;
aqeK.YLabel.Position(1) = -110;
aqeK.YLabel.Position(2) = -3;

% Subfigure labels 
text(axg, -97, 1500, '(a)', 'FontSize',14,'Interpreter','latex','BackgroundColor','white');
text(aqeK, -97, 5,'(b)', 'FontSize',14,'Interpreter','latex','BackgroundColor','white');
text(aqe, -97, 30, '(c)', 'FontSize',14,'Interpreter','latex','BackgroundColor','white');
aqe.Legend.FontSize=16;

% Space-time plot of q_E
figure(3); clf; set(gcf(),"Position", [1 1 800 600])
axt = init_axes(subplot(3,1,[2 3]),'$x$ / km','$t$ / kyr',Ld_km*[0 xmax]);
grid(axt,'minor');
axt.YLabel.Position(1) = -170;

colormap(axt,cmocean('balance'));
qE_max = 50;
% Cut off high and low values for visual distinction
qE_st(qE_st>qE_max) = qE_max; qE_st(qE_st<-qE_max) = -qE_max;

% Refine the mesh grids to reduce aliasing patterns
xmesh = interp2(xmesh, 1);
tmesh = interp2(tmesh, 1);
qE_st = interp2(qE_st,1);

pcolor(axt, Ld_km*xmesh', Td_kyr*tmesh', qEd_mmyr*qE_st');
axt.Children.EdgeColor = 'none';

% Plot grounding line
plot(axt,Ld_km*[xg(1) xg(kplots)] , [-100 Td_kyr*t(kplots)], 'k-');
scatter(axt, Ld_km*xg(kplots), Td_kyr*t(kplots), 32, 'k', 'filled');

% Colour bar
cb = colorbar(axt, 'SouthOutside','FontSize',14,'TickLabelInterpreter','latex');
cb.Label.String = '$q_E$ / mm yr$^{-1}$';
cb.Label.FontSize = 16;
cb.Label.Interpreter = 'latex';

% Separate plot of basin geometry
abg = init_axes(subplot(3, 1, 1),'','$z$ / m',Ld_km*[0 xmax]); 
grid(abg, 'minor');

% dH/dx
HH = S(xmax*X) - b(xmax*X);
Hx = (HH(3:end) - HH(1:end-2)) / (2*xmax*dX*Ld_km);

Hxmax = 9/Ld_km;
Hx(Hx>Hxmax) = Hxmax;
Hx(Hx<-Hxmax) = -Hxmax;

% Colour the basin according to -dH/dx
[m1, m2] = meshgrid(xmax*X(2:end-1), [min(b(xmax*X))-0.1 max(S(xmax*X))+0.1]);
pcolor(abg, Ld_km*m1, Hd*m2, -(Hx'*[1 1])');
abg.Children(1).EdgeColor = 'none';
colormap(abg, cmocean('balance'));

% Draw the basin on
plot(abg, Ld_km*xmax*X, Hd*b(xmax*X), 'k-','linewidth', 1.5);
plot(abg, Ld_km*xmax*X, Hd*S(xmax*X), 'k-','linewidth', 1.5);

text(abg, 50, -4000, '(a)', 'FontSize',14,'EdgeColor','none','BackgroundColor','white',...
    'HorizontalAlignment','Center','Fontsize',14,'interpreter','latex');
text(axt, 50, -10, '(b)', 'FontSize',14,'EdgeColor','none','BackgroundColor','white',...
    'HorizontalAlignment','Center','Fontsize',14,'interpreter','latex');

% Colour bar for -dH/dx
cb2 = colorbar(abg, 'NorthOutside','FontSize',14,'TickLabelInterpreter','latex');
cb2.Label.String = '$- \partial H / \partial x$';
cb2.Label.FontSize = 16;
cb2.Label.Interpreter = 'latex';

