function [R, stable, firstCol] = vdg_routh(c)
%VDG_ROUTH  Routh array of the polynomial with coefficients c = [a_n ... a_0].
%
%   [R, stable, firstCol] = VDG_ROUTH(c)
%
%   R        full Routh array (rows s^n ... s^0, zero padded)
%   stable   true if all first-column entries are strictly positive
%   firstCol first column of the array

c = c(:).';
n = numel(c) - 1;              % polynomial degree

% split into even / odd power coefficients, highest power first
ev = c(1:2:end);
od = c(2:2:end);
if numel(ev) < numel(od), ev(end+1) = 0; end
if numel(od) < numel(ev), od(end+1) = 0; end
ncol = numel(ev);

R = zeros(n+1, ncol);
R(1,:) = ev;
R(2,1:numel(od)) = od;

for i = 3:n+1
    a = R(i-2,:);
    b = R(i-1,:);
    for j = 1:ncol-1
        num = b(1)*a(j+1) - a(1)*b(j+1);
        den = b(1);
        if abs(den) < eps
            den = eps;                 % avoid division by zero
        end
        R(i,j) = num/den;
    end
end

firstCol = R(:,1);
stable = all(firstCol > 0) && all(abs(firstCol) > 0);
end
