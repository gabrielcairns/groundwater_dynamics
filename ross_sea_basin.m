% Ross sea sedimentary basin
% Supplementary animations S6-11
options = optimoptions('fsolve','Display','off');

if ~isfolder('ross_init'); mkdir('ross_init'); end
if ~isfolder('ross_results'); mkdir('ross_results'); end

list = {'1e-10', '1e-12', '3e-13', '1e-13', '3e-14','3e-12'};
[indx,tf] = listdlg('ListString',list);
if tf ==0; error("No experiment selected"); end

gifLength = 15; % in seconds
max_cycles = 20; % how many cycles are allowed in one run, to avoid infinite looping

switch indx
    case 1 
        k0 = 1e-10;  % permeability in m^2
        J = 100;  % x steps
        tres1 = 200; % Time resolution for slow stage
        tres2 = 500; % Time resolution for fast stage
        gifname = 'S7.gif';
    case 2 
        k0 = 1e-12;
        J = 150; 
        tres1 = 300;
        tres2 = 1000; 
        gifname = 'S9.gif';
    case 3
        k0 = 3e-13; 
        J = 100;
        tres1 = 300;
        tres2 = 1000; 
        gifname = 'S10.gif';
    case 4
        k0 = 1e-13; 
        J = 100;  
        tres1 = 300;
        tres2 = 1000; 
        gifname = 'S11.gif';
    case 5 
        k0 = 3e-14;  
        J = 100;
        tres1 = 300;
        tres2 = 1000; 
        gifname = 'S12.gif';
     case 6 
        J = 150;
        k0 = 3e-12; 
        tres1 = 300;
        tres2 = 1000; 
        gifname = 'S8.gif';
end
ICfilename = sprintf('ross_init/ross_init_%g.mat',k0);

alph = 5;
[delta_rho, z_sea, Ld_km, Hd, Td_kyr, r_i, k0_to_Kay, qEd_mmyr] = default_params();
Kay = k0_to_Kay*k0;

% Meshgrid 
X0 = 0;     XJ = 1;
dX = (XJ - X0)/J;       X = (X0:dX:XJ)' ; 
XH = 0.5*(X(1:J)+X(2:J+1));

% Load the basement geometry data
load('bHdata.mat','St_avg','bt_avg','xi_transect');
xgmin = 150/Ld_km;   xgmax = 1640/Ld_km;  xgc = 550/Ld_km;
xslw = 425 /Ld_km;

b = @(x) interp1(xi_transect/Ld_km, bt_avg/Hd, x); 
S = @(x) interp1(xi_transect/Ld_km, St_avg/Hd, x);
H = @(x) S(x) -b(x);

K = 0.8*tres1 + 0.2*tres2;
period = 100/Td_kyr;  % Dimensionless period of cycles

% Create the time mesh
t0 = -period;   tK = 0;     
%t = [reshape((-NCT:-1)+[0:(1/tspc1):(0.8-1/tspc1) 0.8:(1/tspc2):(1-1/tspc2)]',[1,K]) tK];
t = [0:(1/tres1):(0.8-1/tres1) 0.8:(1/tres2):1]-1;
dt = t(2:end) - t(1:end-1);
percents = (K/10):(K/10):K;

% Grounding line (a piecewise function of time)
tc = t/period ;
xg = (((xgmax-xgc)/.8)*(tc+1)+xgc).*(tc<-.2) +...
    (((xgmin-xgmax)/.1)*(tc+.2)+xgmax).*(tc>=-.2).*(tc<-.1) +...
    (((xgc-xgmin)/.1)*(tc+.1)+xgmin).*(tc>=-.1); 
xg_t = (((xgmax-xgc)/.8)*(tc<-.2) +...
    ((xgmin-xgmax)/.1)*(tc>=-.2).*(tc<-.1) +...
    ((xgc-xgmin)/.1)*(tc>=-.1))./period; 

kgif =  [1:1:(0.2*K+1) flip(K+1:-2:(0.2*K+2))];
dtgif = t(kgif(2:end))-t(kgif(1:end-1));
kplots = round(interp1(t,1:K+1,[0.4 0.8 0.85 0.9 0.975 1]-1));

% Ice sheet
PS = @(X,xg) r_i*SSA_ice(X*xg', S(X*xg'), alph);  
   
xmax = max(xi_transect)/Ld_km; 
hplot = zeros(J+1, size(kplots,2));
hgif = zeros(J+1, size(kgif,2));
    
%% 2. The simulation itself
% If no initial condition has been provided using a previous solution, use
% initial saturation with saltwater
total_cycles = 0;
periodic_found = false;
while ~periodic_found && total_cycles < max_cycles
    if ~isfile(ICfilename)
        h_curr = H(xg(1)*X);
    else 
        load(ICfilename,'h_curr');
    end
    h_init = h_curr;        
    u = zeros(J,1);    v = 0.1*ones(J,1);
    exf = zeros(1,K+1); % Contains exit flags of fsolve
    
    qE = zeros(J-1,length(hgif));
    % Values at SLW site
    sslw = S(xslw)*ones(1,length(kgif));
    qEslw = zeros(1,length(kgif));
    
    % Capture intitial value
    if any((1==kplots));    hplot(:, 1==kplots)=h_curr;   end
    if any((1==kgif));      hgif(:, 1==kgif)=h_curr;      end
    for k=1:K
        huv_in = [h_curr;  u;  v]; % Starting guess is previous solution

        ps = PS(X,xg(k+1)); 
        bb = b(xg(k+1)*X);      SS = S(xg(k+1)*X);
        pix = (ps(2:J+1)-ps(1:J)+SS(2:J+1)-SS(1:J))/dX;
     
        f = @(HH) evolve_h_func(HH, h_curr, ps, bb, SS, xg(k+1), xg_t(k), dX, J, delta_rho, Kay, dt(k));
        [huv_out,~,ex]= fsolve(f,huv_in,options);
        
        % Progress updates
        exf(k)=ex;
        if ex<1;    warning('Convergence failed at step %g',k);      end
        if any(k==percents);    fprintf('%g%% complete \n',100*k/K);     end
       
        % Update h
        h_curr = huv_out(1:J+1);    
    
        % Capture h and q_E
        if any((k+1==kplots))    
            hplot(:, k+1==kplots) = h_curr; 
        end
        if any((k+1==kgif))    
            hgif(:, k+1==kgif) = h_curr;     
            G = H(xg(k+1)*X(1:J)).*pix/xg(k+1) + ...
                0.5*(h_curr(2:J+1)+h_curr(1:J)).*(delta_rho*(bb(2:J+1)-bb(1:J)+h_curr(2:J+1)-h_curr(1:J))/(dX*xg(k+1)));
            qE(:,k+1==kgif) = Kay * (G(2:J)-G(1:J-1))/(dX*xg(k+1));
    
            % Values at SLW site
            if xg(k+1)>xslw
                sslw(k+1==kgif) = b(xslw)+interp1(xg(k+1)*X,hgif(:,k+1==kgif),xslw);
                qEslw(k+1==kgif) = interp1(xg(k+1)*X(2:J),qE(:,k+1==kgif),xslw);
            end
        end
        u = huv_out(J+2:2*J+1);    v = huv_out(2*J+2:3*J+1); 
    end
    
    if ~any(exf<1) || find(exf<1,1,'last')<kgif(1)
        disp('Convergence at all stages in cycle')
    else
        warning('Some steps failed to converge in cycle')
    end   
       
    h_final = h_curr;
    % Save current solution as initial condition for future iterations
    save(ICfilename,'h_curr');
    total_cycles = total_cycles + 1;
    periodic_found = (norm(h_init-h_final) < J*1e-4);
    if ~periodic_found
        disp('Solution not periodic; simulating another cycle')
    end
end
if ~periodic_found
    warning("Maximum no. of cycles exceeded, but a periodic solution was not found")
end
% Save the results once out of loop
save(sprintf('ross_results/ross_%g.mat',k0),'Kay','hplot','kgif','xg','J','X','h_curr','qE','t','dX','kplots','PS')

%% 3. gif creation
f1 = figure(2); clf;  set(f1,'Position',[1 1 1000 500],'Color','w');
Movh = init_axes(gca,'$x$ / km', '$z$ / m',Ld_km*[0 xgmax],Hd*[min(b(xmax*X)) max(S(xmax*X)+PS(X,xgmax)/r_i)+0.1]);
plt(length(kgif)) = struct('cdata', [], 'colormap', []);

% Time between sucessive frames
delayT = dtgif*gifLength;
delayT = [delayT delayT(end)];
if min(delayT)<0.02;    delayT =  delayT*0.02/min(delayT);  end
% Make sure framerate is not too fast

% Create each frame
for k=kgif
    cla(Movh);
    title(Movh, sprintf('$\\alpha = %g, k = %g$ m$^2, K =%.3g$, $t=%.0f$ kyr before present',...
        alph,k0,Kay,-Td_kyr*t(k)),'fontsize',14,'Interpreter','latex');
    
    [pice, pint] = draw_frame(Movh,X,b,H,hgif(:, k==kgif),PS(X,xg(k)),xg(k),xgmax);

    write_frame(f1,k,min(kgif),gifname,delayT(k==kgif))
end