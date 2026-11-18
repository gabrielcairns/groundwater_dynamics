% Plotting and processing Antarctic data
% Creates aquifer geometry data along transect
% Fig. 7
% Uses cmocean colormaps
% Uses the following packages:
addpath('./Antarctic Mapping Tools', './Bedmap2 Toolbox for Matlab', './Antarctic Boundaries');
% Greene, C. A., et al. (2016) "Antarctic Mapping Tools for Matlab."
% Fretwell, P., et al. (2013) "Bedmap2: improved ice bed, surface and thickness datasets for Antarctica"
% Rignot, E., et al. (2013) "Ice-shelf melting around Antarctica"
% Also uses basement geometry data from
% Tankersley, D., et al. (2022) "Basement topography and sediment thickness beneath Antarctica's Ross Ice Shelf"
% And grounding line minimum data, reverse engineered by digitising a plot in
% Kingslake, J., et al. (2018) "Extensive retreat and re-advance of the West Antarctic Ice Sheet during the Holocene"

% Import base data (see Tankersley et al. (2022) above)
base = 'Ross_Embayment_basement_filt.nc' ;
% sed = 'Ross_Embayment_sediment.nc' ;

% Subregion of grid
xmax = 5e2; xmin = -8e2;
ymax = -2e2; ymin = -2.2e3;

% Read data (rescale northing/easting to km rather than m)
x = ncread(base, 'x')/1000;       y = ncread(base, 'y')/1000;
[ym, xm] = meshgrid(y((y<ymax)&(y>ymin)), x((x<xmax)&(x>xmin)));
% b
bxy = ncread(base, 'z');     
% Hxy = ncread(sed, 'z'); 
% The Tankersley et al. dataset includes S and H if desired, which in turn
% are based on BedMachine data:
% Morlighem, M., et al. (2020) "Deep glacial troughs and stabilizing ridges 
% unveiled beneath the margins of the Antarctic ice sheet"
bxy = bxy((x<xmax)&(x>xmin), (y<ymax)&(y>ymin));

% S (see Fretwell et al. (2013) above)
[xbm, ybm, Sxy] = bedmap2_data('bed', 1000*xm, 1000*ym, 'xy');
Sxy = interp2(xbm, ybm, Sxy, 1000*xm, 1000*ym);
Hxy = Sxy - bxy;
% Sxy = bxy + Hxy;
Hxy(Hxy<0) = 0; 
% This is not physical, but arises because of the differing datasets

[x_gl, y_gl] = antbounds_data('gl', 'xy');
x_gl = x_gl/1000;     y_gl = y_gl/1000;
[x_shelf_E, y_shelf_E] = antbounds_data('Ross East', 'xy');
x_shelf_E = x_shelf_E/1000;   y_shelf_E = y_shelf_E/1000;
[x_shelf_W, y_shelf_W] = antbounds_data('Ross West', 'xy');
x_shelf_W = x_shelf_W/1000;   y_shelf_W = y_shelf_W/1000;

% Transect (6-point linear spline)
x_points = [-5.6 -4.3 -3 -.7 .4 -.2 -1.2];
y_points = [-3.2 -5.6 -5.5 -6.5 -8.6 -13.3 -16.2];
x_transect = 1e2*x_points;     y_transect = 1e2*y_points;
x_transect = interp1(x_transect, 1:0.01:7);
y_transect = interp1(y_transect, 1:0.01:7);

% Distance along transect; S and b
xi_transect = cumsum([0 sqrt( (x_transect(2:end)-x_transect(1:end-1)).^2 + (y_transect(2:end)-y_transect(1:end-1)).^2)]);
S_transcect = interp2(xm', ym', Sxy', x_transect, y_transect);
b_transect = interp2(xm', ym', bxy', x_transect, y_transect);

% Extrapolating b in the unknown region
firstb = find(~isnan(b_transect), 1);
b_transect(1) = S_transcect(1);
b_transect(isnan(b_transect)) = xi_transect(isnan(b_transect))*(b_transect(firstb)-b_transect(1))/(xi_transect(firstb)) +b_transect(1);

% Perturbations about transect
xtp = zeros(length(x_transect), 5); 
ytp = zeros(length(y_transect), 5);
xitp = zeros(length(y_transect), 5);
Stp = zeros(length(x_transect), 5); 
btp = zeros(length(x_transect), 5);

% Data along the perturbed transects
prtbs = [0 1 2 -1 -2];
for k = 1:5
    % Sampling along perturbed transect
    xtps = 1e2*( x_points + prtbs(k)*[0 0.1 0 0.1 0.1 0.1 0.1]);
    ytps = 1e2*( y_points + prtbs(k)*[0 0.1 0.1 0.1 0 -0.1 -0.1]);    
    xtp(:, k) = interp1(xtps, 1:0.01:7);
    ytp(:, k) = interp1(ytps, 1:0.01:7);
    xitp(:, k) = cumsum([0; sqrt( (xtp(2:end, k)-xtp(1:end-1, k)).^2 + (ytp(2:end, k)-ytp(1:end-1, k)).^2)]);
    Stp(:, k) = interp2(xm', ym', Sxy', xtp(:, k), ytp(:, k)); 
    btp(:, k) = interp2(xm', ym', bxy', xtp(:, k), ytp(:, k));
    
    % Extrapolate btra to the end
    firstb = find(~isnan(btp(:, k)), 1);
    btp(1, k) = Stp(1, k);
    btp(isnan(btp(:, k)), k) = xitp(isnan(btp(:, k)), k)*(btp(firstb, k)-btp(1, k))/(xitp(firstb, k)) +btp(1, k);
    Stp(:, k) = interp1(xitp(:, k), Stp(:, k), xi_transect);
    btp(:, k) = interp1(xitp(:, k), btp(:, k), xi_transect);
end

% Transect averages
St_avg = .2*sum(Stp, 2); bt_avg = .2*sum(btp, 2);

%% Creating plots
% Create figure and axes
figure(1);  clf; set(gcf(),'Position', [ 1 1 1250 575])
ab = init_axes(subplot("Position",[0.42 0.1 0.21 0.83]),'Easting / km','Northing / km',[xmin xmax],[ymin ymax]); axis equal;
aH = init_axes(subplot("Position",[0.74 0.1 0.21 0.83]),'Easting / km','Northing / km',[xmin xmax],[ymin ymax]); axis equal;
aS = init_axes(subplot("Position",[0.1 0.1 0.21 0.83]),'Easting / km','Northing / km',[xmin xmax],[ymin ymax]); axis equal;
% Label fontsize 16

colormap(aS, cmocean('ice'));
colormap(aH, cmocean('matter'));
colormap(ab, flip(cmocean('deep')));

axcell = {ab,aH,aS};

% Plot data
contourf(aS, xm', ym', Sxy', -1600:100:0, 'LineColor', 'none');
contourf(aH, xm', ym', Hxy', 0:1000:8000, 'LineColor', 'none');
contourf(ab, xm', ym', bxy', 'LineColor', 'none');

% Colorbars
cbS = colorbar(aS, 'NorthOutside', 'FontSize',14,'TickLabelInterpreter','latex');
cbS.Label.String = '$S$ / m';
cbS.Label.FontSize = 16;
cbS.Label.Interpreter = 'latex';
cbH = colorbar(aH, 'NorthOutside', 'FontSize',14,'TickLabelInterpreter','latex');
cbH.Label.String = '$H$ / m';
cbH.Label.FontSize = 16;
cbH.Label.Interpreter = 'latex';
cbH.Limits(1)=0;
cbb = colorbar(ab, 'NorthOutside', 'FontSize',14,'TickLabelInterpreter','latex');
cbb.Label.String = '$b$ / m';
cbb.Label.FontSize = 16;
cbb.Label.Interpreter = 'latex';


% Minimum grounding line
% (see Kingslake et al. (2018) above)
kgl = load('kingslake_gl_min.csv');

for ind = 1:3
    % Edge of continental shelf
    contour(axcell{ind}, xm', ym', Sxy', [-1600 -1600], 'k:', 'linewidth', 1);
    % Ice shelf
    plot(axcell{ind}, x_shelf_E, y_shelf_E, 'w:', 'linewidth', 1);
    plot(axcell{ind}, x_shelf_W, y_shelf_W, 'w:', 'linewidth', 1);
    % Current and minimum grounding line
    pglcurr = plot(axcell{ind}, x_gl, y_gl, 'k-', 'linewidth', 1);
    pglmin = plot(axcell{ind}, kgl(:, 1), kgl(:, 2), 'k--', 'linewidth', 2);
    % Onset of sedimentary basin 
    scatter(axcell{ind}, -5.65e2, -3.26e2, 64, 'filled', 'r');
    % Transect
    ts = plot(axcell{ind}, xtp, ytp, 'k-.', 'linewidth', 0.8); ts = ts(1);
    % SLW site
    SLW = plot(aS, -2.8e2, -5.6e2, 'rx', 'markersize', 18, 'linewidth', 3);

    xlim(axcell{ind}, [xmin xmax]);
    ylim(axcell{ind}, [ymin ymax]);
end

% Second figure / axis for plotting transect
figure(2); clf;
atra = init_axes(gca(),'Distance along transect / km', '$z$ / m');

% Colours
cS = cmocean('ice',17);
cb = flip(cmocean('deep',9));
cbt = cb(8,:);
cba = cb(6,:);
cSt = cS(10,:);
cSa = cS(4,:);

plot(atra, xi_transect, Stp, 'Color', cSt, 'linewidth', 2.5);
plot(atra, xi_transect, btp, 'Color', cbt, 'linewidth', 2.5);

% S and b
sav = plot(atra, xi_transect, St_avg, 'Color', cSa, 'linewidth', 2.5);
bav = plot(atra, xi_transect, bt_avg, 'Color', cba, 'linewidth', 2.5);

% Positions of key points
xig_min = interp1(x_transect(1:300), xi_transect(1:300), -490);
xi_SLW = interp1(x_transect(1:300), xi_transect(1:300), -280);
xig_c = interp1(x_transect(1:300), xi_transect(1:300), -170);
xi_end = 1630;

% Plot key points
Ons = scatter(atra, 0, bt_avg(1), 64, 'filled', 'r');
plot(atra, xig_min*[1 1], [-6000 0], 'k--', 'linewidth', 1);
plot(atra, xi_SLW, -640, 'rx', 'markersize', 24, 'linewidth', 3);
plot(atra, xig_c*[1 1], [-6000 0], 'k-', 'linewidth', .75);
pcs = plot(atra, xi_end*[1 1], [-6000 0], 'k:', 'linewidth', 1);

% Legend
legend(atra, [Ons, ts, sav, bav, pglmin, pglcurr, pcs, SLW], ...
    {'Onset of basin', 'Transect', '$S$ (average)', '$b$ (average)', ...
    'Minimum GL', 'Present GL', 'Edge of continental shelf', 'MT measurement'}, ...
    'location', 'southwest', 'interpreter', 'latex', 'fontsize', 14, 'NumColumns',2);

% Subplot labels
text(aS, -750, -300, '(a)', 'FontSize', 16, 'Interpreter', 'Latex', 'BackgroundColor', 'w');
text(ab, -750, -300, '(b)', 'FontSize', 16, 'Interpreter', 'Latex', 'BackgroundColor', 'w');
text(aH, -750, -300, '(c)', 'FontSize', 16, 'Interpreter', 'Latex', 'BackgroundColor', 'w');
text(atra, 20, -800, '(d)', 'FontSize', 16, 'Interpreter', 'Latex', 'BackgroundColor', 'w');

%% Save the final data
save('bHdata.mat', 'xi_transect', 'bt_avg', 'St_avg');
