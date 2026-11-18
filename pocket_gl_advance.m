% Grounding line advance, different K, in parallel
% Figure 6
% Uses cmocean colormap
options = optimoptions('fsolve','Display','off');

[delta_rho, z_sea, Ld_km, Hd, Td_kyr, r_i, k0_to_Kay, qEd_mmyr] = default_params();

Kay = [10 0.5];
L = length(Kay);

% Meshgrid
J = 150;     X0 = 0;     XJ = 1;
dX = (XJ - X0)/J;       X = (X0:dX:XJ)' ; 

% Basement geometry
b = @(x) -.75 -1.5+exp(-(x-.25).^2/.075^2);
H = @(x) -.75 - b(x);

% Create the time mesh
K = 1200;    % # of timesteps
kgif = 1:4:K+1; % subset used for gif
gifLength = 10; % in seconds
t0 = 0;   tK = 8;        
dt = (tK - t0)/K;       t = t0:dt:tK ;
tkyr = Td_kyr*t;
kplots = find(ismember(t,[0 1 5 8])); % where we plot s
% Progress / numerics logging
percents = (K/5):(K/5):K; % for progress logging
exf = zeros(1,K); % Contains exit flags of fsolve

% Grounding line 
xg = .2 + .2*min(t,5);  xgt = 0 + .2*(t<=5);  
xmax = max(xg);

% Ice sheet
PS = @(X,xg) r_i*SSA_ice(X*xg', H(X*xg')+b(X*xg'), .1);

  
%% 2. The solution itself
% Basic steady state (advanced GL)
x = X*xg(end);
ps = PS(X,xg(end));
SS = (ps(end)-ps-H(x)+(1+delta_rho)*(H(x(end))+b(x(end))-b(x)) )/delta_rho;
SS(SS<0)=0;

% Largest pocket steady state (advanced GL)
takeOff = ps(2:end) + H(x(2:end))+(1+delta_rho)*b(x(2:end))>...
    ps(1:end-1) + H(x(1:end-1))+(1+delta_rho)*b(x(1:end-1));
xp = x(find(takeOff,1,'last'));
if isempty(xp); xp=0; end
SSalt = SS;
SSalt(x<xp) = (ps(x==xp)-ps(x<xp)+H(xp)-H(x(x<xp))...
    +(1+delta_rho)*(b(xp)-b(x(x<xp))) )/delta_rho;
SSalt(SSalt<0)=0;
        
% Steady state (initial GL)
x_init = X*xg(1); ps_init = PS(X,xg(1));
SS_init = (ps_init(end)-ps_init-H(x_init)+(1+delta_rho)*(H(x_init(end))+b(x_init(end))-b(x_init)) )/delta_rho;
SS_init(SS_init<0)=0;
% Salt volume compared to final gl position

h_curr = SS_init*ones(1,L); % start in SS
Vf = zeros(L,K+1); % Freshwater volume
hplot = zeros(J+1,L,length(kplots)); % h for plots
hgif = zeros(J+1,L,length(kgif)); % h for making .gif

u = ones(J,2);    v = zeros(J,2);

% Record intial h
Vf(:,1) = sum(H(xg(1)*X) - h_curr,1)*dX*xg(1);
if any((1==kplots));  hplot(:, :,1==kplots) = h_curr; end
if any((1==kgif));  hgif(:, :,1==kgif) = h_curr; end

huv_out = zeros(3*J+1,2);
for k=1:K
    % pS, b, S at each step
    ps = PS(X,xg(k+1)); 
    bb = b(xg(k+1)*X);  S = bb + H(xg(k+1)*X);
    pix = (ps(2:J+1)-ps(1:J)+S(2:J+1)-S(1:J))/dX;

    huv_in = [h_curr;  u;  v]; 
    exf(k) = 3;
    for l = 1:L
        f = @(HH) evolve_h_func(HH, h_curr(:,l), ps, bb, S, xg(k+1), xgt(k), dX, J, delta_rho, Kay(l), dt);
        [huv_out(:,l),~,ex]= fsolve(f,huv_in(:,l),options);
        h_curr(:,l) = huv_out(1:J+1,l);
        exf(k) = min(ex,exf(k));
    end

    % Update progress
        if ex<1;    warning('Convergence failed at step %g',k);      end
        if any(k==percents);    fprintf('%g%% complete \n',100*k/K);     end
    
    % Record data
    Vf(:,k+1) = sum(H(xg(k+1)*X) - h_curr,1)*dX*xg(k+1);
    if any((k+1==kplots))  
        hplot(:, :, k+1==kplots) = h_curr; 
    end
    if any((k+1==kgif))
        hgif(:, :, k+1==kgif) = h_curr;   
    end
end

%% 3. Plotting the results

xtra = [xmax*X; 1.3];
xmaxx = max(xtra);

% Create figure / axes
figure(1); clf; set(gcf(),'Position',[1 1 1000 650])
axg = init_axes(subplot(3,2,5),'$t$ / kyr','$x_g$ / km',Td_kyr*t(kgif([1 end])));     grid(axg,'on'); 
avf = init_axes(subplot(3,2,6),'$t$ / kyr','$\Delta V_f$ / km$^2$',Td_kyr*t(kgif([1 end])));      grid(avf,'on');

letters = 'abcd';
% Plot the grounding line, with snapshot markers
plot(axg, Td_kyr*t(kgif), Ld_km*xg(kgif), 'k-', 'linewidth', 2);
scatter(axg, Td_kyr*t(kplots), Ld_km*xg(kplots), 50, 'black', 'filled');
for P=1:4
    text(axg,Td_kyr*(t(kplots(P))), Ld_km*(xg(kplots(P))-0.18),['(' letters(P) ')'],...
        'EdgeColor','black','BackgroundColor','white','FontSize',14,'interpreter','latex');
end
% Plot Delta Vf
plot(avf, Td_kyr*t(kgif), Ld_km*(Vf(1,kgif)-Vf(2,kgif)), 'r-', 'LineWidth',2);

% Labels and limits
text(axg,10,500,'(e)','fontsize',14,'interpreter','latex','backgroundColor','w');
text(avf,10,100,'(f)','fontsize',14,'interpreter','latex','backgroundColor','w');

% Plot snapshots of s(x,t) at points in time
clrs = cmocean('matter',6);
for P = 1:4
    ah = init_axes(subplot(3,2,P),'$x$ / km','$z$ / m',Ld_km*[0 xmaxx], Hd*[min(b(xmax*X)) max(H(xmax*X)+b(xmax*X)+PS(X,xmax))+0.2]);
    
    % Limits and labels
    title(ah, sprintf('(%s) $t=%g$ kyr',letters(P),round(Td_kyr*t(kplots(P)))),'FontSize',14,'interpreter','latex');
    
    xgkp = xg(kplots(P));
    % Basin, ice, ocean
    plot(ah, Ld_km*[flip(xtra); xtra], Hd*[b(1.3); flip(b(xmax*X)); H(xmax*X)+b(xmax*X); H(1.3)+b(1.3)], 'k-' , 'linewidth', 2);  
    plot(ah, Ld_km*xgkp*X, Hd*(H(xgkp*X)+b(xgkp*X)+ PS(X,xgkp)/r_i), 'c-','linewidth',2);     
    plot(ah, Ld_km*(xgkp+[0 0.07 0.07 0]), ...
        Hd*(H(xgkp)+b(xgkp)+[0 0.2 0.95 1]*(1+delta_rho)*(-H(xgkp)-b(xgkp))/r_i), 'c-','linewidth',2);
    plot(ah, Ld_km*[xgkp+0.07 xmaxx], [0 0], 'b-','linewidth',2);
 
    % s(x,t) for each K and steady state
    K1 = plot(ah, Ld_km*xgkp*X, Hd*(hplot(:,1,P)+b(xgkp*X)), ...
        'color', clrs(2,:), 'linestyle', '-','linewidth',2.5);
    K2 = plot(ah, Ld_km*xgkp*X, Hd*(hplot(:,2,P)+b(xgkp*X)), ...
        'color', clrs(5,:), 'linestyle', '-.','linewidth',2);
    ss = plot(ah, Ld_km*xg(end)*X, Hd*(SSalt+b(xg(end)*X)), 'k--','linewidth',1);
    % Legend
    if P==1 
        legend(ah,[K2 K1 ss],...
            {sprintf('$K=%g$',Kay(2)),sprintf('$K=%g$',Kay(1)),...
            'Final maximal S.S.'},'FontSize',14,'interpreter','latex','location','northeast'); 
    end

end

%% 4. gif creation
% Create figure and axes
gifname = 'S6.gif';
f1= figure(2);  clf;     set(gcf(),'Position',[1 1 900 750],'Color','w');
Movh1 = init_axes(subplot(2,1,1),'$x$ / km', '$z$ / m', Ld_km*[0 xmaxx], Hd*[min(b(xmax*X))-.2 max(H(xmax*X)+b(xmax*X)+PS(X,xmax)/r_i)+0.1]);
Movh2 = init_axes(subplot(2,1,2),'$x$ / km', '$z$ / m', Ld_km*[0 xmaxx], Hd*[min(b(xmax*X))-.2 max(H(xmax*X)+b(xmax*X)+PS(X,xmax)/r_i)+0.1]);

plt(length(kgif)) = struct('cdata', [], 'colormap', []);
delayT = max(gifLength / length(kgif),0.02);

% Create each frame
for k = kgif
    cla(Movh1);     cla(Movh2);
    title(Movh1, sprintf('$\\alpha=0.1, K=%g, t=%g$ kyr before present',Kay(1),800-round(Td_kyr*t(k))),'FontSize',14,'interpreter','latex');
    title(Movh2, sprintf('$\\alpha=0.1, K=%g, t=%g$ kyr before present',Kay(2),800-round(Td_kyr*t(k))),'FontSize',14,'interpreter','latex');
     
    [pice1, pint1] = draw_frame(Movh1,X,b,H,hgif(:, 1,k==kgif),PS(X,xg(k)),xg(k),xmaxx);
    [pice2, pint2] = draw_frame(Movh2,X,b,H,hgif(:, 2,k==kgif),PS(X,xg(k)),xg(k),xmaxx);
    
    write_frame(f1,k,min(kgif),gifname,delayT)
end

