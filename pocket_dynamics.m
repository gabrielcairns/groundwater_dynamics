% Experiment 1: Long-time simulation of interface relaxing to (and past) a pocket steady
% state, using variable timesteps
% Experiment 2: Test the stability, perturbing from a pocket free state
% Fig. 5a and %b, Supplementary animation S4 and S5
% Uses cmocean colormap
options = optimoptions('fsolve','Display','off');

list = {'Long simulation','Stability test'};
[indx,tf] = listdlg('ListString',list);
if tf ==0; error("No experiment selected"); end

Kay = 10;
[delta_rho, z_sea, Ld_km, Hd, Td_kyr, r_i, k0_to_Kay, qEd_mmyr] = default_params();

% Meshgrid
if indx == 1 
    J = 150;  
elseif indx==2 
    J = 200; 
end   
X0 = 0;     XJ = 1;
dX = (XJ - X0)/J;       X = (X0:dX:XJ)' ; 

% Topography
b = @(x) -2.5+exp(-(x-.25).^2/.025^2);
H = @(x) -1 - b(x);

% Create the time mesh
if indx == 1
    t = [0:0.001:.1 .104:0.004:.3 .31:.01:1 1.04:.04:3 3.1:.1:10 10.4:.4:50 51:1:100];
    dt = t(2:end)-t(1:end-1);
    K = length(t)-1;    % # of timesteps
    % When to record
    kplots = find(ismember(t,[.1 1 10 100])); % plotting stage
    kgif = 1:2:K+1;
    gifLength = 6; % in seconds
elseif indx == 2
    t = 0:0.005:1;
    dt = t(2:end)-t(1:end-1);
    K = length(t)-1;    % # of timesteps
    % Timesteps at which to capture h(x,t)
    kplots = find(ismember(t,[0 t(end)]));
    kgif = 1:1:K+1;
    gifLength = 5; % in seconds
end
percents = (K/5):(K/5):K; % Progress updates

% Grounding line
xg = .5*ones(size(t));  xgt = zeros(size(t));
xmax = max(xg); 

% Ice sheet
%PS = @(X,xg) 0r_i*viscousIce(X*xg', H(X*xg')+b(X*xg'), .05);
PS = @(X,xg) r_i*SSA_ice(X*xg', H(X*xg')+b(X*xg'), .05);

hplot = zeros(J+1,size(kplots,2));
hgif = zeros(J+1,size(kgif,2));
Vleft = zeros(1,size(kplots,2)); % Volume to the left of x_p^*    

%% 2. The simulation itself
% Find steady state
xe = X*xg(end); 
ps = PS(X,xg(end));
SS = (ps(end)-ps-H(xe)+(1+delta_rho)*(H(xe(end))+b(xe(end))-b(xe)) )/delta_rho;
xn = xe(find(SS<=0,1,'last'));
SS(SS<0)=0;

% Find alternative steady state (largest pocket)
takeOff = ps(2:end) + H(xe(2:end))+(1+delta_rho)*b(xe(2:end))>...
    ps(1:end-1) + H(xe(1:end-1))+(1+delta_rho)*b(xe(1:end-1));
xet = xe(takeOff);
xp = xet(end);
SS_alt = SS;
SS_alt(xe<xp) = (ps(xe==xp)-ps(xe<xp)+H(xp)-H(xe(xe<xp))...
    +(1+delta_rho)*(b(xp)-b(xe(xe<xp))) )/delta_rho;
SS_alt(SS_alt<0) = 0;

% Initial condition: h=H
if indx == 1
    h_curr = H(xg(1)*X);
elseif indx == 2
    h_curr = SS + 0.1*tanh((xg(1)-xe)/.05);
end

u = ones(J,1);    v = zeros(J,1);
exf = zeros(1,K); % Contains exit flags of fsolve

% Record initial condition
if any((1==kplots))  
    hplot(:, 1==kplots) = h_curr; 
    Vleft(:, 1==kplots) = sum(h_curr.*(xg(1)*X<xp))*xg(1)*dX; 
end
if any((1==kgif))    
    hgif(:, 1==kgif) = h_curr; 
end
for k=1:K
    %fun = FUN(hplus, hI, PS, b, S, XG, XGT, PARAMS)
    hplusI = [h_curr;  u;  v]; % Starting guess is previous solution

    ps = PS(X,xg(k+1)); 
    bb = b(xg(k+1)*X);  S = bb +H(xg(k+1)*X);
    pix = (ps(2:J+1)-ps(1:J)+S(2:J+1)-S(1:J))/dX;
    f = @(HH) evolve_h_func(HH, h_curr, ps, bb, S, xg(k+1), xgt(k), dX, J, delta_rho, Kay, dt(k));

    [huv_out,~,ex]= fsolve(f,hplusI,options);
    
    % Progress updates
    exf(k)=ex;
    if ex<1;    warning('Convergence failed at step %g',k);      end
    if any(k==percents);    fprintf('%g%% complete \n',100*k/K);     end
    
    % Update and record h
    h_curr=huv_out(1:J+1);
    if any((k+1==kplots))  
        hplot(:, k+1==kplots) = h_curr; 
        Vleft(:, k+1==kplots) = sum(h_curr.*(xg(k+1)*X<xp))*xg(k+1)*dX; 
    end
    if any((k+1==kgif))    
        hgif(:, k+1==kgif) = h_curr; 
    end

    % Update complementarity variables
    u = huv_out(J+2:2*J+1);    v = huv_out(2*J+2:3*J+1); 
end

%% 3. Plotting the results
x_plus = 0.6;
xtra = [xmax*X; x_plus];
clrs = cmocean('matter',6);

% Still image
figure(1); clf; set(gcf(),'Position', [1 1 1050 375]);
ah = init_axes(gca(),'$x$ / km','$z$ / m',Ld_km*[0 max(xtra)],Hd*[min(b(xmax*X)) max(H(xmax*X)+b(xmax*X)+PS(X,xmax))+0.1]); 
% Plot outline of sedimentary basin
plot(ah, Ld_km*[flip(xtra); xtra], Hd*[b(0.6); flip(b(xmax*X)); H(xmax*X)+b(xmax*X); H(0.6)+b(0.6)], 'k-' , 'linewidth', 2);  

if indx == 1
% Plot h at different times
pice = plot(ah, Ld_km*xg(k)*X, Hd*(hplot(:,1)+b(xg(kplots(1))*X)), 'Color', clrs(2,:), 'LineStyle', '-.', 'LineWidth', 2);
pint = plot(ah, Ld_km*xg(k)*X, Hd*(hplot(:,2)+b(xg(kplots(2))*X)), 'Color', clrs(3,:), 'LineStyle', '-', 'LineWidth', 2);
pSS = plot(ah, Ld_km*xg(k)*X, Hd*(hplot(:,3)+b(xg(kplots(3))*X)), 'Color', clrs(4,:), 'LineStyle', '-.', 'LineWidth', 2);
p4 = plot(ah, Ld_km*xg(k)*X, Hd*(hplot(:,4)+b(xg(kplots(4))*X)), 'Color', clrs(5,:), 'LineStyle', '-', 'LineWidth', 2);

% Plot steady state
pSS=  plot(ah, Ld_km*xg(end)*X, Hd*(SS_alt+b(xg(end)*X)), 'k--','linewidth',1);
% Ice sheet
pp = plot(ah, Ld_km*xmax*X, ...
     Hd*(H(xmax*X)+b(xmax*X)+PS(X,xmax)/r_i), 'c-','linewidth',2);
plot(ah, Ld_km*(xmax+[0 0.07 0.07 0]), ...
    Hd*(H(xmax)+b(xmax)+[0 0.2 0.95 1]*(41/40)*(-H(xmax)-b(xmax))/r_i), 'c-','linewidth',2.5);
plot(ah, Ld_km*[0.57 x_plus], [0 0], 'b-','linewidth',2.5);

% Legend/ labels
legend(ah,[pice,pint,pSS,p4,pSS,pp],...
    {'$t=10^1$ kyr', '$t=10^2$ kyr', '$t=10^3$ kyr', '$t=10^4$ kyr', 'Maximal pocket', 'Ice surface'},...
    'fontsize',14,'interpreter','latex','location','southeast');
text(ah, 5, -1200, '(a)', 'fontsize',14,'interpreter','latex');

elseif indx == 2
    fill(ah, Ld_km*[xet; flip(xet)], Hd*[H(xet)+b(xet); b(flip(xet))], [0.75 0.75 0.75], 'EdgeColor','none');

% Plot h at different times
pstart = plot(ah, Ld_km*xg(1)*X, Hd*(hplot(:,1)+b(xg(kplots(1))*X)),...
    'Color', clrs(2,:), 'linestyle', '-.', 'linewidth', 3);
pend = plot(ah, Ld_km*xg(3)*X, Hd*(hplot(:,2)+b(xg(kplots(2))*X)),...
    'Color', clrs(5,:), 'linestyle', '--', 'linewidth', 3);
% Plot pocket steady state
pSS = plot(ah, Ld_km*xg(end)*X, Hd*(SS_alt+b(xg(end)*X)), 'k--','linewidth',1.5);
% x_p
scatter(ah, Ld_km*xet([1 end]), Hd*b(xet([1 end])), 60, 'filled', 'k');

% Plot p at different times
pp = plot(ah, Ld_km*xmax*X, Hd*(H(xmax*X)+b(xmax*X)+PS(X,xmax)/r_i), 'c-','linewidth',2.5);
plot(ah, Ld_km*(xmax+[0 0.07 0.07 0]), ...
    Hd*(H(xmax)+b(xmax)+[0 0.2 0.95 1]*(41/40)*(-H(xmax)-b(xmax))/r_i), 'c-','linewidth',2.5);
plot(ah, Ld_km*[0.57 x_plus], [0 0],'b-','linewidth',2.5);

% Legend/ labels
legend(ah,[pstart,pend,pSS,pp],...
    {'$t=0$ kyr',sprintf('$t=%g$ kyr',Td_kyr*t(end)), 'Maximal pocket', 'Ice surface'},...
    'FontSize', 14, 'interpreter', 'latex', 'location', 'southeast');
text(ah,5,-1200,'(b)','fontsize',14,'interpreter','latex');

% Label the regions with arrows
text(ah, 23, -1200, '$\mathcal{U}>0$', 'fontsize', 14, 'interpreter', 'latex')
text(ah, 103, -1200, '$\mathcal{U}<0$', 'fontsize', 14, 'interpreter', 'latex')
text(ah, 143, -1200, '$\mathcal{U}>0$', 'fontsize', 14, 'interpreter', 'latex')
text(ah, 25, -1400, '$\rightarrow$', 'fontsize', 28, 'interpreter', 'latex')
text(ah, 105, -1400, '$\leftarrow$', 'fontsize', 28, 'interpreter', 'latex')
text(ah, 145, -1400, '$\rightarrow$', 'fontsize', 28, 'interpreter', 'latex')
text(ah, Ld_km*xp, Hd*(b(xp)-.15), '$x_p^*$', 'fontsize', 14, 'interpreter', 'latex')
end

%% 4. gif creation
f1 = figure(2); clf;  set(gcf(),'Position',[1 1 900 350],'Color','w');
Movh = init_axes(gca,'$x$ / km', '$z$ / m',Ld_km*[0 x_plus],Hd*[min(b(xmax*X))-.2 max(H(xmax*X)+b(xmax*X)+PS(X,xmax))+0.1]);
if indx == 1 
    gifname = 'S4.gif'; 
elseif indx == 2 
    gifname = 'S5.gif'; 
end

plt(length(kgif)) = struct('cdata', [], 'colormap', []);
delayT = max(gifLength / length(kgif),0.02);

% Create each frame
for k=kgif
    cla(Movh);
    if indx ==1
        title(Movh, sprintf('$\\alpha=0.05, K=%g, t=%.0f$ kyr',Kay,Td_kyr*(100-t(k))), 'Interpreter', 'latex');
    elseif indx == 2
         title(Movh, sprintf('$\\alpha=0.05, K=%g, t=%.0f$ kyr before present',Kay,Td_kyr*(1-t(k))), 'fontsize', 14, 'Interpreter', 'latex');
    end

    [pice, pint] = draw_frame(Movh,X,b,H,hgif(:, k==kgif),PS(X,xg(k)),xg(k),x_plus);

    if indx ==2
        pSS=plot(Movh, Ld_km*xg(end)*X, Hd*(SS_alt+b(xg(end)*X)), 'k--','linewidth',2);
        legend(Movh,[pice,pint,pSS],{'$S(x)+H_i(x,t)$','$s(x,t)$','Steady state'},'fontsize',14,'interpreter','latex','location','southeast');
    end
    
    write_frame(f1,k,min(kgif),gifname,delayT)
end