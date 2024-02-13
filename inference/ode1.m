%% FUNCTION Y = ODE1(PARS)
%  Input model parameters, pars, output a vector y that corresponds to the
%  "experimental" data:
%
%   Y = [c(1),c(2),...,c(28),i(1),i(2),...,i(28)]
%
%  where c(t) is the cell count under continuous treatment at time t
%  [days], and i(t) is that for cells under intermittent treatment.
%
%  For this model, pars corresponds to:
%
%   pars = [g1,g2,g3,g4,omega]
%
%  For simplicity, we assume that the initial condition is measured exactly
%  (given as an optional parameter).
function Y = ode1(pars,varargin)

    % Get initial condition
    if nargin == 2
        n0 = varargin{1};
    else
        n0 = 1;
    end

    % Get parameters
    g1 = pars(1); g2 = pars(2); g3 = pars(3); g4 = pars(4);
    omega = pars(5);

    % Growth rate function (effectively the continuous PFM equivalent)
    gr_drugon  = @(x) g2 + x * (g4 - g2);
    gr_drugoff = @(x) g1 + x * (g3 - g1);
    gr = @(x,drug) gr_drugon(x) * drug + gr_drugoff(x) * ~drug;

    % Time points (exclude zero in output)
    T = 1:1:28;

    % Drug on or off
    drug_cont = @(t) true;
    drug_int  = @(t) mod(t,14) < 7;

    % Model
    function dy = rhs(t,y,drugfcn)
        n = y(1); x = y(2); drug = drugfcn(t);

        dn = n * gr(x,drug);
        dx = omega * drug - omega * ~drug;
        if drug && x >= 1
            dx = 0;
        elseif ~drug && x <= 0
            dx = 0;
        end
        dy = [dn;dx];

    end

    % Determine stop points (i.e., so we don't overshoot x = 0 or x = 1 by too much)
    tst = sort([T, (0:7:28) + 1 / omega + 1e-10]);

    % Solve (continuous therapy)
    [~,Ycont] = ode45(@(t,y) rhs(t,y,drug_cont),tst,[n0;0]);

    % Solve (intermittent therapy)
    [~,Yint]  = ode45(@(t,y) rhs(t,y,drug_int),tst,[n0;0]);

    % Return
    Y = [interp1(tst,Ycont(:,1),T),interp1(tst,Yint(:,1),T)];

end