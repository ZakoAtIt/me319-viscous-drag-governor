function vdg_main()
%VDG_MAIN  Reproduce every number and figure reported in the paper.
%
%   Run from the sim/ directory:      vdg_main
%   or from the repository root:      run('sim/vdg_main.m')
%
%   Produces figures in ../figs/ and a machine-readable summary in
%   ../docs/results.txt.  Uses only base MATLAB / GNU Octave functions.

close all; clc;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
figdir  = fullfile(here, '..', 'figs');
docdir  = fullfile(here, '..', 'docs');
if ~exist(figdir,'dir'), mkdir(figdir); end
if ~exist(docdir,'dir'), mkdir(docdir); end

p = vdg_params();
mp = vdg_model(p);            % proportional design
mi = vdg_model(p, p.Ti);      % proportional + reset (PI) design

fid = fopen(fullfile(docdir,'results.txt'),'w');
logf = @(varargin) both(fid, varargin{:});

logf('==================================================================\n');
logf(' ME319 Assignment 1 - Viscous-Drag Governor (VDG)\n');
logf(' Reproduction of all reported results.  MATLAB/Octave: %s\n', version);
logf('==================================================================\n\n');

%% ---------------- 1. Plant / controller / loop transfer functions ------
logf('--- Transfer functions (deviation variables) ---\n');
logf('Gp(s) = %s / %s\n', poly2str(mp.Gpnum), poly2str(mp.Gpden));
logf('Gc(s) = %s / %s   [P design]\n', poly2str(mp.Cnum), poly2str(mp.Cden));
logf('Gc(s) = %s / %s   [PI design]\n', poly2str(mi.Cnum), poly2str(mi.Cden));
logf('Loop L(s) numerator   (P) = %s\n', poly2str(mp.Lnum));
logf('Loop L(s) denominator (P) = %s\n\n', poly2str(mp.Lden));

%% ---------------- 2. Closed-loop poles and stability ------------------
logf('--- Characteristic polynomial and closed-loop poles ---\n');
logf('P  design: char(s) = %s\n', poly2str(mp.char));
rp = roots(mp.char);
logf('P  poles: %s\n', mat2str(rp,6));
[R, st] = vdg_routh(mp.char);
logf('P  Routh first column: %s\n', mat2str(R(:,1).',6));
logf('P  Routh-Hurwitz stable: %d\n', st);

logf('PI design: char(s) = %s\n', poly2str(mi.char));
ri = roots(mi.char);
logf('PI poles: %s\n', mat2str(ri,6));
[Ri, sti] = vdg_routh(mi.char);
logf('PI Routh first column: %s\n', mat2str(Ri(:,1).',6));
logf('PI Routh-Hurwitz stable: %d\n\n', sti);

% dominant pair of the P design = the roots closest to the imaginary axis
[~, islow] = max(real(rp));
dom = rp(islow);
wn = abs(dom); zeta = -real(dom)/wn;
logf('P dominant pole %.6f  ->  wn = %.4f rad/s, zeta = %.4f\n', dom, wn, zeta);
if imag(dom) ~= 0
    logf('  (closed-form overshoot of an equivalent 2nd-order pair: %.2f %%)\n', ...
         100*exp(-pi*zeta/sqrt(1-zeta^2)));
end
logf('P droop (steady-state speed error per unit load): %.4f (%.2f %%)\n', ...
     mp.droop, 100*mp.droop);
logf('\n');

%% ---------------- 3. Frequency response margins ------------------------
w = logspace(-3,2,40000);
[PMp, wcp, GMp, wgp] = margins(mp, w);
[PMi, wci, GMi, wgi] = margins(mi, w);
logf('--- Stability margins ---\n');
logf('P : gain-crossover %.4f rad/s,  PM = %.2f deg,  GM = %.2f dB @ %.4f rad/s\n', ...
     wcp, PMp, 20*log10(GMp), wgp);
logf('PI: gain-crossover %.4f rad/s,  PM = %.2f deg,  GM = %.2f dB @ %.4f rad/s\n\n', ...
     wci, PMi, 20*log10(GMi), wgi);

%% ---------------- 4. Time-domain simulation ----------------------------
T = 80; dt = 2e-3;
dP = 0.30;   % +30 % load step
dM = -0.30;  % -30 % load step

% open loop:  J dw/dt = -b dw - dTL,  throttle held at xv0
Aol = -p.b/p.J; Bol = -1/p.J;
[t, XolP] = vdg_rk4(Aol, Bol, 1, dP, T, dt);
[~, XolM] = vdg_rk4(Aol, Bol, 1, dM, T, dt);
wolP = p.w0 + XolP(:,1); wolM = p.w0 + XolM(:,1);
xvOL = ones(size(t));

% closed loop, P and PI
[~, XcP, YcP] = vdg_rk4(mp.A, mp.B, mp.C, dP, T, dt);
[~, XcM, YcM] = vdg_rk4(mp.A, mp.B, mp.C, dM, T, dt);
[~, XiP, YiP] = vdg_rk4(mi.A, mi.B, mi.C, dP, T, dt);
[~, XiM, YiM] = vdg_rk4(mi.A, mi.B, mi.C, dM, T, dt);

res = struct();
res.olP = metrics(t, wolP, xvOL, 1);
res.olM = metrics(t, wolM, xvOL, 1);
res.clP = metrics(t, p.w0 + YcP(:,1), p.xv0 + YcP(:,2), p.w0);
res.clM = metrics(t, p.w0 + YcM(:,1), p.xv0 + YcM(:,2), p.w0);
res.piP = metrics(t, p.w0 + YiP(:,1), p.xv0 + YiP(:,2), p.w0);
res.piM = metrics(t, p.w0 + YiM(:,1), p.xv0 + YiM(:,2), p.w0);

logf('--- Transient response metrics ---\n');
logf('%-30s %9s %10s %9s %9s %9s %9s\n', 'case','w_inf','ss_err[%]','peak[%]','ts2%[s]','ts.5%[s]','xv_inf');
dump(logf, '+30% load step', res);
dump(logf, '-30% load step', res);
logf('\n');

%% ---------------- 5. Figures -------------------------------------------
% One colour per case, everywhere:  blue = VDG proportional,  green = P + reset,
% grey dashed = open loop.  The speed of the disturbance is carried by the panel
% title, never by the colour.
blue  = [0 0.24 0.60];
green = [0 0.50 0.15];
grey  = [0.40 0.40 0.40];

% ---- Fig 3: open vs closed loop, with a dedicated zoom panel -----------
% Colour identity is GLOBAL and never changes panel to panel:
%   grey dashed = open loop,  blue = VDG closed loop.
% The SIGN of the load step is carried by the panel title (wide panels) or by
% the LINE STYLE (solid = +30%, dash-dot = -30%), never by a change of colour.
% Every panel carries its own legend, placed OUTSIDE the axes, so that no
% legend box ever covers a curve.
lg = {'Location','northoutside','Orientation','horizontal','FontSize',6.5,'Box','off'};

f = figure('Position',[100 100 900 330],'Color','w');
set(f,'DefaultAxesFontSize',9);
kz = t <= 6;  tk = t(kz);
wzP = p.w0 + YcP(kz,1);  wzM = p.w0 + YcM(kz,1);
posS = {[0.055 0.40 0.255 0.52],[0.370 0.40 0.255 0.52],[0.685 0.40 0.295 0.52]};
posX = [0.055 0.075 0.925 0.235];

% (row 1, col 1) speed, +30%
axes('Position',posS{1}); hold on; grid on; box on;
plot(t, wolP,        '--', 'Color', grey, 'LineWidth', 1.2, 'DisplayName','open loop');
plot(t, p.w0+YcP(:,1), '-', 'Color', blue, 'LineWidth', 1.5, 'DisplayName','VDG closed');
plot([0 T],[1 1], 'k:', 'HandleVisibility','off');
plot([0 T],[res.clP.w_inf res.clP.w_inf], ':', 'Color', blue, 'HandleVisibility','off');
text(30, res.clP.w_inf-0.055, sprintf('\\omega_\\infty=%.4f', res.clP.w_inf), 'FontSize',7);
xlim([0 T]); ylim([0.66 1.06]); xlabel('t [s]'); ylabel('\omega(t)');
title('+30% load step'); legend(lg{:});

% (row 1, col 2) speed, -30%
axes('Position',posS{2}); hold on; grid on; box on;
plot(t, wolM,        '--', 'Color', grey, 'LineWidth', 1.2, 'DisplayName','open loop');
plot(t, p.w0+YcM(:,1), '-', 'Color', blue, 'LineWidth', 1.5, 'DisplayName','VDG closed');
plot([0 T],[1 1], 'k:', 'HandleVisibility','off');
plot([0 T],[res.clM.w_inf res.clM.w_inf], ':', 'Color', blue, 'HandleVisibility','off');
text(30, res.clM.w_inf+0.045, sprintf('\\omega_\\infty=%.4f', res.clM.w_inf), 'FontSize',7);
xlim([0 T]); ylim([0.94 1.34]); xlabel('t [s]'); ylabel('\omega(t)');
title('-30% load step (disengaged)'); legend(lg{:});

% (row 1, col 3) speed, zoom: the peak and the settling time are visible here
axes('Position',posS{3}); hold on; grid on; box on;
plot(tk, 1.02*ones(size(tk)), 'k--', 'DisplayName','\pm2% band');
plot(tk, 0.98*ones(size(tk)), 'k--', 'HandleVisibility','off');
plot(tk, wzP, '-',  'Color', blue, 'LineWidth', 1.6, 'DisplayName','+30% (solid)');
plot(tk, wzM, '-.', 'Color', blue, 'LineWidth', 1.6, 'DisplayName','-30% (dash-dot)');
plot([res.clP.ts res.clP.ts],[0.975 1.025], ':', 'Color', grey, 'HandleVisibility','off');
[ypk, ipk] = max(abs(wzP - 1));
plot(tk(ipk), wzP(ipk), 'o', 'Color', blue, 'MarkerFaceColor','w', 'HandleVisibility','off');
text(tk(ipk)+0.30, wzP(ipk)-0.013, sprintf('peak %.2f%%', 100*ypk), 'FontSize',7);
text(res.clP.ts+0.20, 0.9762, sprintf('t_s=%.2f s', res.clP.ts), 'FontSize',7);
xlim([0 6]); ylim([0.975 1.025]); xlabel('t [s]'); ylabel('\omega(t)');
title('zoom on the governed transient'); legend(lg{:});

% (row 2) throttle valve: ONE wide panel, so the trace stays legible
axes('Position',posX); hold on; grid on; box on;
plot(t, p.xv0+YcP(:,2), '-',  'Color', blue, 'LineWidth', 1.5, 'DisplayName','VDG, +30% (solid)');
plot(t, p.xv0+YcM(:,2), '-.', 'Color', blue, 'LineWidth', 1.5, 'DisplayName','VDG, -30% (dash-dot)');
plot(t, xvOL,           '--', 'Color', grey, 'LineWidth', 1.2, 'DisplayName','open loop (throttle held)');
plot([0 T],[1 1], 'k:', 'HandleVisibility','off');
xlim([0 T]); ylim([0.63 1.36]); xlabel('t [s]'); ylabel('x_v(t)');
legend(lg{:});
export(f, figdir, 'fig_open_closed');

% ---- Fig: zoom on the closed-loop transient ---------------------------
f = figure('Position',[100 100 560 360],'Color','w'); hold on; grid on; box on;
k = t <= 10;
plot(t(k), p.w0+YcP(k,1), '-', 'Color', blue,  'LineWidth', 1.8, 'DisplayName','+30% load step');
plot(t(k), p.w0+YcM(k,1), '-', 'Color', green, 'LineWidth', 1.8, 'DisplayName','-30% load step');
plot([0 10],[1.02 1.02], 'k--', 'HandleVisibility','off');
plot([0 10],[0.98 0.98], 'k--', 'HandleVisibility','off');
plot([0 10],[res.clP.w_inf res.clP.w_inf], ':', 'Color', blue, 'HandleVisibility','off');
xlim([0 10]); ylim([0.965 1.035]);
xlabel('t [s]'); ylabel('\omega(t)');
legend('Location','northeast','FontSize',8);
title('Closed-loop speed regulation (VDG, proportional)');
export(f, figdir, 'fig_zoom');

% ---- Fig 4: proportional vs reset (PI) --------------------------------
% NOTE: this figure is included at \columnwidth, so its on-screen size is kept
% small (roughly column width) to stop the labels shrinking to illegibility.
f = figure('Position',[100 100 490 235],'Color','w');
set(f,'DefaultAxesFontSize',9);
subplot(1,2,1); hold on; grid on; box on;
plot(t, p.w0+YcP(:,1), '-', 'Color', blue,  'LineWidth', 1.5, 'DisplayName','P (droop)');
plot(t, p.w0+YiP(:,1), '-', 'Color', green, 'LineWidth', 1.5, 'DisplayName','P + reset');
plot([0 40],[1 1], 'k:', 'HandleVisibility','off');
plot([0 40],[1.005 1.005], 'k:', 'HandleVisibility','off');
xlim([0 40]); ylim([0.970 1.010]); xlabel('t [s]'); ylabel('\omega(t)');
title('+30% load step');
legend(lg{:});
subplot(1,2,2); hold on; grid on; box on;
plot(t, p.w0+YcM(:,1), '-', 'Color', blue,  'LineWidth', 1.5, 'DisplayName','P (droop)');
plot(t, p.w0+YiM(:,1), '-', 'Color', green, 'LineWidth', 1.5, 'DisplayName','P + reset');
plot([0 40],[1 1], 'k:', 'HandleVisibility','off');
plot([0 40],[0.995 0.995], 'k:', 'HandleVisibility','off');
xlim([0 40]); ylim([0.990 1.030]); xlabel('t [s]'); ylabel('\omega(t)');
title('-30% load step (machine disengaged)');
legend(lg{:});
export(f, figdir, 'fig_reset');

% ---- Fig: pole map ----------------------------------------------------
f = figure('Position',[100 100 640 380],'Color','w'); hold on; grid on; box on;
plot(real(rp), imag(rp), 'x', 'Color', blue,  'MarkerSize', 10, 'LineWidth', 1.8, 'DisplayName','VDG (P)');
plot(real(ri), imag(ri), 'o', 'Color', green, 'MarkerSize', 8,  'LineWidth', 1.6, 'DisplayName','VDG (P + reset)');
plot([0 0],[-10 10],'k-', 'HandleVisibility','off');
xlabel('Re(s) [rad/s]'); ylabel('Im(s) [rad/s]');
legend('Location','southwest','FontSize',8);
title('Closed-loop pole locations'); xlim([-10 1]); ylim([-6 6]);
export(f, figdir, 'fig_poles');

% ---- Fig: loop Bode with margins --------------------------------------
Lp = freqresp(mp, w); Li = freqresp(mi, w);
f = figure('Position',[100 100 620 520],'Color','w');
set(f,'DefaultAxesFontSize',9);
subplot(2,1,1); hold on; grid on; box on;
semilogx(w, 20*log10(abs(Lp)), 'Color', blue,  'LineWidth', 1.6, 'DisplayName','P');
semilogx(w, 20*log10(abs(Li)), 'Color', green, 'LineWidth', 1.4, 'DisplayName','P + reset');
plot([1e-3 10],[0 0], 'k--', 'HandleVisibility','off');
xlim([1e-3 10]); ylabel('|L(j\omega)| [dB]');
legend('Location','southwest','FontSize',8);
subplot(2,1,2); hold on; grid on; box on;
semilogx(w, unwrapdeg(angle(Lp)*180/pi), 'Color', blue,  'LineWidth', 1.6);
semilogx(w, unwrapdeg(angle(Li)*180/pi), 'Color', green, 'LineWidth', 1.4);
plot([1e-3 10],[-180 -180], 'k--', 'HandleVisibility','off');
xlim([1e-3 10]); xlabel('\omega [rad/s]'); ylabel('\angle L(j\omega) [deg]');
export(f, figdir, 'fig_bode');

logf('Figures written to: %s\n', figdir);
fclose(fid);
fprintf('\nDone. Results: %s\n', fullfile(docdir,'results.txt'));
type(fullfile(docdir,'results.txt'));
end

%% ======================= helper functions ==============================
function m = metrics(t, wv, xv, ref)
m.w_inf  = wv(end);
m.ss_err = 100*abs(wv(end)-ref)/ref;
m.peak   = 100*max(abs(wv-ref))/ref;
m.xv_inf = xv(end);
dfin = abs(wv(end)-ref);
if dfin > 1e-9
    band = 0.02*dfin;
    idx = find(abs(abs(wv-ref) - dfin) > band);
    if isempty(idx), m.ts = 0; else, m.ts = t(min(idx(end)+1, numel(t))); end
else
    m.ts = NaN;
end
idx = find(abs(wv-ref) > 0.005*ref);
if isempty(idx)
    m.ts5 = 0;
elseif idx(end) >= numel(t)
    m.ts5 = NaN;          % never enters the +-0.5% band (permanent offset)
else
    m.ts5 = t(idx(end)+1);
end
m.xv_dev = max(abs(xv-1));
end

function dump(logf, tag, r)
% pick the open-loop case that matches the sign of the disturbance
if contains(tag,'+')
    row(logf, [tag ' open loop'], r.olP, 1.0);
    row(logf, [tag ' VDG (P)'],      r.clP, r.clP.xv_inf);
    row(logf, [tag ' VDG (P+reset)'], r.piP, r.piP.xv_inf);
else
    row(logf, [tag ' open loop'], r.olM, 1.0);
    row(logf, [tag ' VDG (P)'],      r.clM, r.clM.xv_inf);
    row(logf, [tag ' VDG (P+reset)'], r.piM, r.piM.xv_inf);
end
end

function row(logf, tag, m, xv)
logf('%-30s %9.4f %10.3f %9.3f %9s %9s %9.4f\n', tag, m.w_inf, m.ss_err, m.peak, ...
     num2str(m.ts,'%.2f'), num2str(m.ts5,'%.2f'), xv);
end

function [PM, wc, GM, wg] = margins(m, w)
L = freqresp(m, w);
mag = abs(L);
[~, i] = min(abs(mag-1));
wc = w(i); PM = 180 + angle(L(i))*180/pi;
ph = unwrap(angle(L))*180/pi;
j = find(ph <= -180, 1);
if isempty(j)
    GM = Inf; wg = NaN;
else
    GM = 1/mag(j); wg = w(j);
end
end

function L = freqresp(m, w)
s = 1i*w(:).';
L = polyval(m.Lnum, s) ./ polyval(m.Lden, s);
end

function y = unwrapdeg(x)
y = unwrap(x*pi/180)*180/pi;
end

function s = poly2str(c)
c = c(:).';
n = numel(c)-1;
terms = {};
for k = 1:numel(c)
    pw = n-(k-1);
    if abs(c(k)) < 1e-12, continue; end
    if pw == 0
        terms{end+1} = sprintf('%.6g', c(k));                       %#ok<AGROW>
    elseif pw == 1
        terms{end+1} = sprintf('%.6g s', c(k));                     %#ok<AGROW>
    else
        terms{end+1} = sprintf('%.6g s^%d', c(k), pw);              %#ok<AGROW>
    end
end
s = strjoin(terms, ' + ');
end

% exportgraphics() crops tightly to the figure, so the resulting PDF has no
% page-sized white margin when it is scaled into a LaTeX column.  print() is
% kept as a fallback for older MATLAB / Octave.
function export(f, figdir, name)
% Force a LIGHT theme.  In MATLAB R2025a+ the editor theme is inherited by new
% figures, so without this the exported plots come out with a near-black axes
% background and pale grey text, which prints very badly.
set(f,'Color','w');
try
    theme(f,'light');
catch
    % pre-R2025a: set the colours by hand (also keeps Octave happy)
end
ax = findall(f,'Type','axes');
for k = 1:numel(ax)
    a = ax(k);
    try, set(a,'Color','w','XColor','k','YColor','k','ZColor','k'); end
    try, set(a,'GridColor',[0.80 0.80 0.80],'MinorGridColor',[0.90 0.90 0.90]); end
    try, set(a.Title,'Color','k'); end
    try, set(a.XLabel,'Color','k'); set(a.YLabel,'Color','k'); end
    try
        lg = get(a,'Legend');
        if ~isempty(lg)
            set(lg,'Color','w','TextColor','k','EdgeColor',[0.30 0.30 0.30]);
        end
    catch
    end
end
set(f,'Toolbar','none','MenuBar','none');
% Hide the per-axes toolbar as well: exportgraphics otherwise captures the small
% hovering toolbar in the corner of the exported figure.
try
    for k = 1:numel(findall(f,'Type','axes'))
        ax = findall(f,'Type','axes');
        set(ax(k).Toolbar,'Visible','off');
    end
catch
    % older MATLAB releases / Octave have no axes toolbar
end
pdf = fullfile(figdir,[name '.pdf']);  png = fullfile(figdir,[name '.png']);
try
    exportgraphics(f, pdf, 'ContentType','vector', 'BackgroundColor','white');
    exportgraphics(f, png, 'Resolution',300, 'BackgroundColor','white');
catch
    set(f,'PaperPositionMode','auto');
    print(f, '-dpdf', pdf);
    print(f, '-dpng', '-r300', png);
end
end

function both(fid, varargin)
fprintf(varargin{:});
fprintf(fid, varargin{:});
end
