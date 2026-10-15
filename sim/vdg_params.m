function p = vdg_params()
%VDG_PARAMS  Physical parameters of the Viscous-Drag Governor (VDG).
%
%   All quantities are expressed in the normalised system required by the
%   assignment: the nominal engine speed is omega0 = 1 and the nominal load
%   torque is TL0 = 1.  Time is in seconds.
%
%   Dimensional anchors used to justify the numbers (see paper, Table I):
%     nominal speed      n0  = 60 rpm   -> omega0 = 6.283 rad/s
%     nominal power      P0  = 20 kW    -> nominal torque T0 = 3183 N.m
%   A normalised inertia J corresponds to a physical inertia I = J*T0/omega0.

% ---------------- Plant: steam engine + flywheel + throttle valve ----------
p.J     = 10;    % Flywheel + engine inertia.  I = 10*3183/6.283 = 5.1e3 kg.m^2
                 % (a 2 t flywheel of ~1.6 m radius of gyration).  Gives a
                 % spin-down time constant J/b = 10 s, typical of a mill engine.
p.b     = 1;     % Lumped viscous damping (journal friction + velocity-dependent
                 % drag of the belt-driven line shaft).  Chosen so that the
                 % damping torque at the nominal point equals the nominal load
                 % torque (b*omega0 = TL0 = 1).  This reproduces the classically
                 % documented behaviour of an ungoverned engine, whose speed
                 % change is proportional to the load change (1/b = 1).
%   NOTE the nominal torque balance.  The damping torque is INSIDE TL (see the
%   linearisation), so the shaft equation at the nominal point requires
%       Kt*xv0 = TL(omega0) = TL0 = 1
%   and with xv0 = 1 the throttle gain is Kt = 1.  An earlier revision used
%   Kt = 2, which double-counted the damping as a separate resisting torque.
p.Kt    = 1;     % Throttle-to-torque gain (nominal balance Kt*xv0 = TL0 = 1).
p.tau_v = 0.1;   % Throttle-valve time constant (light, stiff valve, 0.1 s).

% ---------------- Controller: viscous drag-cup speed sensor ----------------
%   The three sensor parameters are quoted after dividing the whole drum equation
%   by the spring stiffness ks, so the numbers that set the DYNAMICS are the
%   ratios cs/ks (sensor gain, s), Js/ks (sensor time constant squared) and
%   Kl*cs/ks (controller DC gain).  These ratios are what the paper reports; the
%   absolute values below are rescaled by 1/50 so that the film transmits a
%   physically sensible torque.  Scaling Js, cs and ks together changes NOTHING in
%   the simulation - it only fixes the physical interpretation.
p.Js    = 0.0004; % Sensing-drum inertia; Js/ks = 0.02 s^2  ->  wn_s = 7.07 rad/s.
p.cs    = 0.004;  % Oil-film shear coefficient; cs/ks = 0.2 s.
                  % At nominal speed the film transmits cs*omega0 = 12.7 N.m, so it
                  % dissipates 12.7*6.283 = 80 W, i.e. 0.4% of the 20 kW engine.
                  % (A 50x larger film would give 637 N.m and 4 kW - 20% of engine
                  % power - which would also break the isothermal assumption.)
p.ks    = 0.02;   % Torsion-spring stiffness; theta_ss = (cs/ks)*omega0 = 0.2
                  % normalised units = 72 deg full scale, and 1.3 deg per 30%
                  % load step - comfortably resolvable.
p.Kl    = 80;     % Linkage gain, dx_v/dtheta.  With Kt = 1 the DC loop gain is
                  % G0 = Kt*Kl*cs/(ks*b) = 16, as in the paper.  Physically a
                  % 0.40 m drum lever feeds a 20:1 compound lever to the valve
                  % stem; a 30% load step moves the stem by 0.28 normalised units
                  % = 28 mm, driven by a 1.3 deg drum rotation.

% ---------------- Optional reset (integral) time constant ------------------
p.Ti    = 4;     % Reset dashpot time constant of the PI variant (s).

% ---------------- Normalised operating point --------------------------------
p.w0  = 1;
p.TL0 = 1;
p.xv0 = 1;

end
