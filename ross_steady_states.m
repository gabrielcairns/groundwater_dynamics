% Plot steady states & potential pockets for a range of alpha
% Fig. 8
% Uses cmocean colormap

alphas = [0.005 0.05 0.5 5];
L = length(alphas);
[delta_rho, z_sea, Ld_km, Hd, Td_kyr, r_i, k0_to_Kay, qEd_mmyr] = default_params();

% Load basement geometry
load('bHdata.mat', 'St_avg', 'bt_avg','xi_transect');

b = @(x) interp1(xi_transect/Ld_km, bt_avg/Hd, x); 
S = @(x) interp1(xi_transect/Ld_km, St_avg/Hd, x);
H = @(x) S(x) -b(x);

% Ice sheet
PS = @(X,xg,alph) r_i*SSA_ice(X*xg', S(X*xg'), alph);

% Meshgrid
dX = 0.001; X = 0:dX:1;
xgs = (150:30:1650)/Ld_km;

k = 49; % Choose advanced grounding line position
XG = xgs(k);
xe = XG*X; 

% Matrix containing ps for each alpha
ps = zeros(L,length(X));
for l = 1:L
    ps(l,:) = PS(X,XG,alphas(l)); 
end

% Basic steady state
SSnw = (ps(:,end)-ps-H(xe)+(1+delta_rho)*(H(xe(end))+b(xe(end))-b(xe)) )/delta_rho;
SSw = SSnw;
SSw(SSw<0)=0;

% Find where steady states are possible
takeOff = ps(:,2:end) + H(xe(2:end))+(1+delta_rho)*b(xe(2:end))>...
    ps(:,1:end-1) + H(xe(1:end-1))+(1+delta_rho)*b(xe(1:end-1));
h_alt = SSw;
% Ignore region where b is extrapolated
takeOff(:,xe<400/Ld_km)=0;
% Working backwards, find potential pockets and include them
for j = length(xe):-1:50
    takeOff(:,(1:length(xe)-1)>=j)=0;
    for l = 1:L
        if any(takeOff(l,:))
            xp = max(xe(takeOff(l,:)));
            h_alt(l,xe<xp) = (ps(l,xe==xp)-ps(l,xe<xp)+H(xp)-H(xe(xe<xp))...
                +(1+delta_rho)*(b(xp)-b(xe(xe<xp))) )/delta_rho;
            h_alt(l,h_alt(l,:)<0)=0;
            takeOff(l,:) = logical(((1:length(xe)-1)<j).*(h_alt(l,2:end)==0).*(ps(l,2:end) + H(xe(2:end))+(1+delta_rho)*b(xe(2:end))>...
                ps(l,1:end-1) + H(xe(1:end-1))+(1+delta_rho)*b(xe(1:end-1))));
            takeOff(l,xe<400/Ld_km)=0;
        end
    end
end

figure(1); clf; set(gcf(),'Position',[1 1 1050 350])
ax = init_axes(gca(),'$x$ / km','$z$ / m',[0 1644],[min(b(xe))*Hd 2300]); 

clrs = cmocean('matter',6);
% Plot each of the steady states
for l = 1:L
    plot(ax,Ld_km*xe, Hd*(b(xe)+h_alt(l,:)), 'Color',clrs(l+1,:),'LineWidth',2);
end

% Plot basement geometry
plot(ax,xi_transect, bt_avg,'k-','LineWidth',2);
plot(ax,xi_transect, St_avg,'k-','LineWidth',2);

% Plot the respective ice sheets
for l = 1:L
   plot(ax,Ld_km*xe, Hd*(S(xe)+ps(l,:)),'Color',clrs(l+1,:),'LineWidth',2);
end

% Ice shelf, ocean
plot(ax,Ld_km*(xe(end))+ 15*[0 1 1 0], Hd*(S(xe(end))+[0 0.2 1 1]*ps(1,end)), 'Color',clrs(5,:),'linewidth',2);    
plot(ax,[Ld_km*(xe(end))+15 1644],[0 0],'k-','LineWidth',2);

% Legend
legend(ax,[ax.Children(3:6)], {'$\alpha=5$','$\alpha=0.5$','$\alpha=0.05$','$\alpha=0.005$'},'fontsize',14,'interpreter','latex','location','southwest');