% Solution for periodic grounding line movement 
% Fig. 3, Supplementary animations S1-3
% Uses cmocean colormap
options = optimoptions('fsolve','Display','off');

if ~isfolder('periodic_init'); mkdir('periodic_init'); end

% Case selection
list = {'high K','med K','low K','very low K'};
[indx,tf] = listdlg('ListString',list);
if tf ==0; error("No experiment selected"); end

[delta_rho, z_sea, Ld_km, Hd, Td_kyr, r_i, k0_to_Kay, qEd_mmyr] = default_params();
alph = .1;

gifLength = 10; % in seconds
max_cycles = 50; % how many cycles are allowed in one run, to avoid infinite looping
switch indx
    case 1
        Kay = 10; 
        init_cond_filename = 'periodic_init/periodic_init_fast.mat';
        letter = '(a)';
        qEYs = [-200 400];
        gifname = 'S1_graphic.gif';
    case 2
        Kay = 1; 
        init_cond_filename = 'periodic_init/periodic_init_med.mat';
        letter = '(b)';
        qEYs = [-25 125];
        gifname = 'S2_graphic.gif';
    case 3
        Kay = 0.1;
        init_cond_filename = 'periodic_init/periodic_init_slow.mat';
        letter = '(c)';
        qEYs = [-25 125];
        gifname = 'S3_graphic.gif';
    case 4
        Kay = 0.01;
        init_cond_filename = 'periodic_init/periodic_init_vslow.mat';
        letter = '';
        qEYs = [-5 40];
        gifname = 'veryslo_graphic.gif';
end
K = 100;

% Meshgrid
J = 100;     X0 = 0;     XJ = 1;
dX = (XJ - X0)/J;       X = (X0:dX:XJ)' ; 

b = @(x) -3.5*ones(size(x));
H = @(x)  2.5*ones(size(x));   
% Stages at which to capture h(x,t)
kplots =  round(linspace(1,K+1,9));
kgif = 1:2:K+1;
ka = kplots(1:5);   kr = kplots(5:9);
% Plots of advance and retreat

% Create the time mesh
t0 = 0; tK = 1;
dt = (tK - t0)/K;       t = t0:dt:tK ;
tkyr = Td_kyr*t;
percents = (K/5):(K/5):K;

xg_max = .3;
xg = xg_max - .1*cos(2*pi*t);
xg_t = .1*(2*pi)*sin(2*pi*t); 

% Ice sheet
PS = @(X,xg) r_i*SSA_ice(X*xg', H(X*xg')+b(X*xg'), alph);

xmax = max(xg);     xmin=min(xg);
hplot = zeros(J+1, size(kplots,2));
hgif = zeros(J+1, size(kgif,2));

%% 2. The simulation itself
% Find steady states (at equal dt timesteps)
SSplot = zeros(J+1,length(ka));
qESS = zeros(J-1,length(ka));
for k = ka
    XG = xg(k);
    xe = X*XG;    ps = PS(X,XG); 
    SSw = (ps(end)-ps-H(xe)+(1+delta_rho)*(H(xe(end))+b(xe(end))-b(xe)) )/delta_rho;
    SSplot(SSw>0,k==ka)=SSw(SSw>0);
    uf = (ps(2:J+1)-ps(1:J) + H(xe(2:J+1))-H(xe(1:J)) + (1+delta_rho)*(b(xe(2:J+1))-b(xe(1:J))))/(XG*dX); 
    G = 0.5*H(xe(2:J+1))+H(xe(1:J)).*uf + delta_rho*0.5*(SSw(2:J+1)+SSw(1:J)).*(SSw(2:J+1)-SSw(1:J))/(XG*dX) ;
    qESS(:,k==ka) = Kay*(G(2:end)-G(1:end-1))/(XG*dX);
end  

total_cycles = 0;
periodic_found = false;
while ~periodic_found && total_cycles < max_cycles
    
    if ~isfile(init_cond_filename)
        h_curr = H(xg(1)*X);
    else 
        load(init_cond_filename,'h_curr');
    end
    
    h_init = h_curr;
    u = ones(J,1);    v = zeros(J,1); % Complementarity variables
    exf = zeros(1,K); % Contains exit flags of fsolve
    qE = zeros(J-1,size(hplot,2));
    
    % Capture initial values
    ps = PS(X,xg(1)); 
    bb = b(xg(1)*X);  S = bb + H(xg(1)*X);
    pix = (ps(2:J+1) - ps(1:J) + S(2:J+1) - S(1:J))/dX;
    if any((1==kplots))
        hplot(:, 1==kplots) = h_curr; 
        uf = pix/xg(1);
        us = uf + delta_rho*(bb(2:J+1)-bb(1:J)+h_curr(2:J+1)-h_curr(1:J))/(dX*xg(1)) ;        
        G = H(xg(1)*X(1:J)).*uf + ...
            0.5*(h_curr(2:J+1)+h_curr(1:J)).*(us-uf);
        qE(:, 1==kplots) = Kay*(G(2:J)-G(1:J-1))/(dX*xg(1));
    end
    if any((1==kgif))
        hgif(:, 1==kgif) = h_curr; 
    end
    
    % Iteration:
    for k=1:K
        huv_in = [h_curr;  u;  v]; 
        % Starting guess is previous solution
        
        ps = PS(X,xg(k+1)); 
        bb = b(xg(k+1)*X);  S = bb + H(xg(k+1)*X);
        pix = (ps(2:J+1) - ps(1:J) + S(2:J+1) - S(1:J))/dX;
    
        f = @(HH) evolve_h_func(HH, h_curr, ps, bb, S, xg(k+1), xg_t(k+1), ...
            dX, J, delta_rho, Kay, dt); % Function to solve
        [huv_out,~,ex] = fsolve(f,huv_in,options);
        
        % Progress updates
        exf(k) = ex;
        if ex<1;    warning('Convergence failed at step %g',k);      end
        if any(k==percents);    fprintf('Cycle %g%% complete \n',100*k/K);     end
        
        % Update h
        h_curr = huv_out(1:J+1);
        
        % Capture current values
        if any((k+1==kplots))  
            hplot(:, k+1==kplots) = h_curr;
            uf = pix/xg(k+1) + delta_rho*(bb(2:J+1)-bb(1:J))/(dX*xg(k+1));
            us = uf + delta_rho*(h_curr(2:J+1)-h_curr(1:J))/(dX*xg(k+1)) ;        
            G = H(xg(k+1)*X(1:J)).*uf + ...
                0.5*(h_curr(2:J+1)+h_curr(1:J)).*(us-uf);
            qE(:,k+1==kplots) = Kay*(G(2:J)-G(1:J-1))/(dX*xg(k+1));
        end
        if any((k+1==kgif))    
            hgif(:, k+1==kgif) = h_curr; 
        end
    
        % Update complementarity variables
        u = huv_out(J+2:2*J+1);    v = huv_out(2*J+2:3*J+1); 
    end

    h_final = h_curr;
    % Save current solution as initial condition for future iterations
    save(init_cond_filename,'h_curr');
    total_cycles = total_cycles + 1;
    periodic_found = (norm(h_init-h_final) < J*1e-4);
    if ~periodic_found
        disp('Solution not periodic; simulating another cycle')
    end
end
if ~periodic_found
    warning("Maximum no. of cycles exceeded, but a periodic solution was not found")
end

%% 3. Plotting the results

figure(1);  clf;     set(gcf(),'Position',[1 1 1050 400]);
Ys = Hd*[min(b(xmax*X)) max(PS(X,xmax)/r_i + H(xmax*X)+b(xmax*X))+0.1];
aha = init_axes(subplot(3,2,[1 3]),'','$z$ / m',Ld_km*[0 0.5],Ys);
ahr = init_axes(subplot(3,2,[2 4]),'','$z$ / m',Ld_km*[0 0.5],Ys);
aqEa = init_axes(subplot(3,2,5),'$x$ / km','$q_E$ / mm yr$^{-1}$', Ld_km*[0 0.5],qEYs);
aqEr = init_axes(subplot(3,2,6),'$x$ / km','$q_E$ / mm yr$^{-1}$',Ld_km*[0 0.5],qEYs);
grid(aqEa,'on');    grid(aqEr,'on');

% Colour map
clrs = cmocean('matter',6);
clrsa =  clrs(2:6,:);
clrsr =  clrs(2:6,:);

% Several plots of h and q_E on one image
plot(aha, Ld_km*xmax*[1.25; flip(X); X; 1.25], ...
    Hd*[b(0.5); flip(b(xmax*X)); H(xmax*X)+b(xmax*X); H(0.5)+b(0.5)], 'k-', 'LineWidth', 2);
plot(ahr, Ld_km*xmax*[1.25; flip(X); X; 1.25], ...
    Hd*[b(0.5); flip(b(xmax*X)); H(xmax*X)+b(xmax*X); H(0.5)+b(0.5)], 'k-', 'LineWidth', 2);
for k = ka % Advance
    plot(aha, Ld_km*xg(k)*X, Hd*(hplot(:,k==kplots)+b(xg(k)*X)),'linestyle' ,'-', 'color', clrsa(k==ka,:), 'LineWidth', 1);
    plot(aqEa, Ld_km*xg(k)*X(2:J), qEd_mmyr*qE(:,k==kplots),'linestyle', '-', 'color', clrsa(k==ka,:), 'LineWidth', 1);
end
pa = aha.Children(3);
for k = kr % Retreat
    plot(ahr, Ld_km*xg(k)*X, Hd*(hplot(:,k==kplots)+b(xg(k)*X)), '-', 'color', clrsr(k==kr,:), 'LineWidth', 1);
    plot(aqEr, Ld_km*xg(k)*X(2:J), qEd_mmyr*qE(:,k==kplots), '-', 'color', clrsr(k==kr,:), 'LineWidth', 1);
end
pr = ahr.Children(3);
% Steady state h
plot(aha, Ld_km*xmax*X, Hd*(SSplot(:,5)+b(xmax*X)), 'k--', 'LineWidth', .5);
plot(aha, Ld_km*xmin*X, Hd*(SSplot(:,1)+b(xmin*X)), 'k--', 'LineWidth', .5);
pSS = plot(aha, Ld_km*xg_max*X, Hd*(SSplot(:,3)+b(xg_max*X)), 'k--', 'LineWidth',  .5);
plot(ahr, Ld_km*xmax*X, Hd*(SSplot(:,5)+b(xmax*X)), 'k--', 'LineWidth', .5);
plot(ahr, Ld_km*xmin*X, Hd*(SSplot(:,1)+b(xmin*X)), 'k--', 'LineWidth', .5);
plot(ahr, Ld_km*xg_max*X, Hd*(SSplot(:,3)+b(xg_max*X)), 'k--', 'LineWidth', .5);

% Steady state q_E
plot(aqEa, Ld_km*xmax*X(2:J), qEd_mmyr*qESS(:,1), 'k--', 'LineWidth', .5);
plot(aqEa, Ld_km*xmin*X(2:J), qEd_mmyr*qESS(:,5), 'k--', 'LineWidth', .5);
plot(aqEa, Ld_km*xg_max*X(2:J), qEd_mmyr*qESS(:,3), 'k--', 'LineWidth', .5);
plot(aqEr, Ld_km*xmax*X(2:J), qEd_mmyr*qESS(:,1), 'k--', 'LineWidth', .5);
plot(aqEr, Ld_km*xmin*X(2:J), qEd_mmyr*qESS(:,5), 'k--', 'LineWidth', .5);
plot(aqEr, Ld_km*xg_max*X(2:J), qEd_mmyr*qESS(:,3), 'k--', 'LineWidth', .5);

% Ice sheet and ocean
plot(aha, Ld_km*xmin*X, Hd*(H(xmin*X)+b(xmin*X)+PS(X,xmin)/r_i), 'c-', 'LineWidth', 2.5);
pp = plot(aha, Ld_km*xmax*X, Hd*(H(xmax*X)+b(xmax*X)+PS(X,xmax)/r_i), 'c-', 'LineWidth', 2.5);
plot(aha, Ld_km*(xmin+[0 0.07 0.07 0]), ...
    Hd*(H(xmin)+b(xmin)+[0 0.2 0.95 1]*(1+delta_rho)*(-H(xmin)-b(xmin))/r_i), 'c-','linewidth',2.5);
plot(aha, Ld_km*(xmax+[0 0.07 0.07 0]), ...
    Hd*(H(xmax)+b(xmax)+[0 0.2 0.95 1]*(1+delta_rho)*(-H(xmax)-b(xmax))/r_i), 'c-','linewidth',2.5);
plot(aha, Ld_km*[xmin+0.07 0.5], [0 0], 'b-', 'LineWidth', 2.5);
plot(ahr, Ld_km*xmin*X, Hd*(H(xmin*X)+b(xmin*X)+PS(X,xmin)/r_i), 'c-', 'LineWidth', 2.5);
pp = plot(ahr, Ld_km*xmax*X, Hd*(H(xmax*X)+b(xmax*X)+PS(X,xmax)/r_i), 'c-', 'LineWidth', 2.5);
plot(ahr, Ld_km*(xmin+[0 0.07 0.07 0]), ...
    Hd*(H(xmin)+b(xmin)+[0 0.2 0.95 1]*(1+delta_rho)*(-H(xmin)-b(xmin))/r_i), 'c-','linewidth',2.5);
plot(ahr, Ld_km*(xmax+[0 0.07 0.07 0]), ...
    Hd*(H(xmax)+b(xmax)+[0 0.2 0.95 1]*(1+delta_rho)*(-H(xmax)-b(xmax))/r_i), 'c-','linewidth',2.5);
plot(ahr, Ld_km*[xmin+0.07 0.5], [0 0], 'b-', 'LineWidth', 2.5);

text(aha, 5, -100, [letter '(i)'], 'fontsize', 14, 'interpreter', 'latex');
scatter(aha, Ld_km*xg_max, Hd*(b(xg_max)+H(xg_max)), 32, 'filled', 'k');
text(aha, Ld_km*xg_max, -800, '$\bar{x}_g$', 'fontsize', 14, 'interpreter', 'latex');

if indx == 1
legend(aha, [pa,pSS,pp], {'$s(x,t)$', 'Steady state', 'Ice surface'}, ...
    'fontsize', 14, 'interpreter', 'latex', 'location' ,'southeast');
end

text(ahr, 5, -100, [letter '(ii)'], 'fontsize', 14, 'interpreter', 'latex');
scatter(ahr, Ld_km*xg_max, Hd*(b(xg_max)+H(xg_max)), 32, 'filled', 'k');
text(ahr, Ld_km*xg_max, -800, '$\bar{x}_g$', 'fontsize', 14, 'interpreter', 'latex');


text(aqEa,5,qEYs(2)-.2*(qEYs(2)-qEYs(1)), [letter '(iii)'], 'fontsize', 14, 'interpreter', 'latex');
text(aqEr,5,qEYs(2)-.2*(qEYs(2)-qEYs(1)), [letter '(iv)'], 'fontsize', 14, 'interpreter', 'latex');

aha.XTickLabels = '';
ahr.XTickLabels = '';

if indx == 1
    title(aha,'Advance','Fontsize',16,'Interpreter','latex')
    title(ahr,'Retreat','Fontsize',16,'Interpreter','latex')
end

%% 4. gif creation
f1 = figure(2);  clf;   set(gcf(),'Position',[1 1 800 450],'Color','w');
xmax = 0.5;
Movh = init_axes(gca,'$x$ / km', '$z$ / m',Ld_km*[0 xmax], Hd*[min(b(xmax*X))-.2 max(PS(X,xmax)/r_i + H(xmax*X)+b(xmax*X))+0.1]);
plt(length(kgif)) = struct('cdata', [], 'colormap', []);
delayT = max(gifLength / length(kgif),0.02);

% Creating each frame
for k = kgif
    cla(Movh);
    title(Movh,sprintf('$K=%g, \\alpha=%g, t=%.0f$ kyr before present',Kay,alph,Td_kyr*(1-t(k))),'FontSize',14,'Interpreter','Latex');
    
    [pice, pint] = draw_frame(Movh,X,b,H,hgif(:, k==kgif),PS(X,xg(k)),xg(k),xmax);
    % Steady states
    pSS=plot(Movh, Ld_km*xg(end)*X, Hd*(SSplot(:,1)+b(xg(end)*X)),'k--', 'linewidth', 2);
    plot(Movh, Ld_km*.4*X, Hd*(SSplot(:,5)+b(.4*X)),'k--', 'linewidth', 2);
    plot(Movh, Ld_km*.3*X, Hd*(SSplot(:,3)+b(.3*X)),'k--', 'linewidth', 2)

    legend(Movh,[pice,pint,pSS],{'$S(x)+H_i(x,t)$','$s(x,t)$','Steady state'},...
        'fontsize',14,'Interpreter','Latex','Location','southeast');    
        
    write_frame(f1,k,min(kgif),gifname,delayT);

end