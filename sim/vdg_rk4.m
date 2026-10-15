function [t, X, Y] = vdg_rk4(A, B, C, dTL, T, dt)
%VDG_RK4  Integrate x' = A x + B*dTL for a step disturbance of size dTL.
%
%   [t, X, Y] = VDG_RK4(A, B, C, dTL, T, dt)
%
%   Uses a classical fourth-order Runge-Kutta scheme so that no Control
%   System Toolbox (or Octave control package) is required.  Works in both
%   MATLAB and GNU Octave.

if nargin < 6 || isempty(dt), dt = 2e-3; end
if nargin < 5 || isempty(T),  T  = 80;   end

N = round(T/dt);
t = (0:N)'*dt;
n = size(A,1);
X = zeros(N+1, n);

f = @(x) A*x + B*dTL;

x = zeros(n,1);
X(1,:) = x.';
for k = 1:N
    k1 = f(x);
    k2 = f(x + (dt/2)*k1);
    k3 = f(x + (dt/2)*k2);
    k4 = f(x + dt*k3);
    x  = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    X(k+1,:) = x.';
end

Y = (C*X.').';      % outputs: columns [dw, dxv]

end
