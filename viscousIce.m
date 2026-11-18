function Hi = viscousIce(x,bed,alph)
    % Generalised viscous ice sheet model
    % Solve d/dx { Hi^N1 * |dHi/dx + dbed/dx|^N2 } = alph
    % subject to rho_i Hi(xg) = rho_s |bed(xg)| and no flux at x=0
    N1 = 4; N2 = 3; % The viscous coefficients

    % Meshgrid
    dx = x(2:end)-x(1:end-1);
    Hi = ones(size(x));

    % First integration
    q1n = (alph*x).^(1/N2);
    % Boundary condition
    Hi(end) = -bed(end)*(41/40)/.917;
    % Second integration
    for k=length(x):-1:2
        Hi(k-1) = Hi(k)  + bed(k)-bed(k-1)+ dx(k-1) *q1n(k)./(Hi(k)).^(N1/N2);
    end
end
