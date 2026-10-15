function m = vdg_model(p, Ti)
%VDG_MODEL  Linear time-invariant model of plant, controller and closed loop.
%
%   m = VDG_MODEL(p)      proportional (spring-restrained) VDG
%   m = VDG_MODEL(p, Ti)  VDG with viscous reset, i.e. PI action
%
%   Deviation variables about the nominal operating point, state
%       x = [dw, dth, dthd, dxv, z]
%   with
%       J  dw/dt              = -b*dw + Kt*dxv - dTL
%       Js d2th/dt2 + cs dth/dt + ks th = cs*(dw + z)
%       dz/dt                 = dw / Ti          (reset / integral action)
%       tau_v ddxv/dt         = -Kl*dth - dxv
%
%   Inputs :  n/a (only the load disturbance dTL is used)
%   Outputs:  y = [dw ; dxv]
%
%   Normalised throttle position is reported as xv = xv0 + dxv.

if nargin < 2 || isempty(Ti)
    n = 4;
    m.hasReset = false;
else
   n = 5;
    m.hasReset = true;
end

A = zeros(n); B = zeros(n,1); C = zeros(2,n);

% --- engine shaft -------------------------------------------------------
A(1,1) = -p.b/p.J;  A(1,4) = p.Kt/p.J;  B(1,1) = -1/p.J;

% --- viscous drag-cup sensing drum -------------------------------------
A(2,3) = 1;
A(3,1) =  p.cs/p.Js;  A(3,2) = -p.ks/p.Js;  A(3,3) = -p.cs/p.Js;

% --- throttle valve -----------------------------------------------------
A(4,2) = -p.Kl/p.tau_v;  A(4,4) = -1/p.tau_v;

% --- viscous reset (integral) path --------------------------------------
if m.hasReset
    A(3,5) = p.cs/p.Js;      % reset angle added to the drum drive
    A(5,1) = 1/p.Ti;         % integrator of the speed error
end

C(1,1) = 1;   % dw
C(2,4) = 1;   % dxv

m.n = n; m.A = A; m.B = B; m.C = C; m.D = zeros(2,1);

% ---------------- controller C(s) = dxv / dw (deviation variables) -------
% P  : C(s) = Kl*cs / (Js s^2 + cs s + ks)
% PI : C(s) = Kl*cs / (Js s^2 + cs s + ks) * (Ti s + 1)/(Ti s)
m.Cnum = p.Kl*p.cs;
m.Cden = [p.Js p.cs p.ks];
if m.hasReset
    m.Cnum = conv(m.Cnum, [Ti 1]);
    m.Cden = conv(m.Cden, [Ti 0]);
end

% ---------------- plant: u -> dw ----------------------------------------
% Gp(s) = Kt / [ (J s + b)(tau_v s + 1) ]      (valve followed by shaft)
m.Gpnum = p.Kt;
m.Gpden = conv([p.J p.b], [p.tau_v 1]);

% ---------------- loop transfer L0(s) = Gp(s) * C0(s), C0 = +Kl cs/(...) ---
% The governor gain is physically negative (rising speed closes the throttle);
% with dxv = -Kl*theta the closed loop obeys 1 - Gp*Gc = 0, i.e. 1 + L0 = 0
% with L0 > 0.  L0 is therefore the loop transfer in the standard negative-
% feedback form used for the Bode / phase-margin assessment.
m.Lnum = conv(m.Gpnum, m.Cnum);
m.Lden = conv(m.Gpden, m.Cden);

% ---------------- closed-loop characteristic polynomial -----------------
% 1 + L0 = 0  ->  Lden + Lnum = 0
m.char = m.Lden + [zeros(1, numel(m.Lden)-numel(m.Lnum)) m.Lnum];

% ---------------- disturbance transfer dw / dTL -------------------------
% dw/dTL = -Gd/(1+L0)  with Gd = 1/(Js+b), giving
%   dw/dTL = -(tau_v s + 1)*Cden / char(s)
numD = -conv([p.tau_v 1], m.Cden);
m.Gd_num = numD;
m.Gd_den = m.char;

% DC gain of the disturbance path (the droop, speed error per unit load)
m.droop = abs(polyval(numD,0)/polyval(m.char,0));

end
