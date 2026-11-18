function [pice,  pint, pshelf, pbasin, pocean] = draw_frame(ax,X,b,H,h,ps,xg,xgmax)

[drho, ~, Ldkm, Hd, ~, ri, ~, ~] = default_params();

pbasin = plot(ax, Ldkm*[flip(xgmax*X); xgmax*X], Hd*[flip(b(xgmax*X)); b(xgmax*X)+H(xgmax*X)], 'k-','linewidth',2); 
% Ocean and ice
pocean = plot(ax, Ldkm*[xg+0.07 xgmax], [0 0], 'b-','linewidth',2);
pice = plot(ax, Ldkm*xg*X, Hd*(b(xg*X)+H(xg*X)+ ps/ri), 'c-','linewidth',2); 
pshelf = plot(ax, Ldkm*(xg+ 0.07*[0 1 1 0]), Hd*(b(xg)+H(xg)+[0 0.2 0.95 1]*(1+drho)*(-b(xg)-H(xg))/ri), 'c-','linewidth',2);
pint = plot(ax, Ldkm*xg*X, Hd*(h+b(xg*X)), 'r-', 'linewidth', 2);

end