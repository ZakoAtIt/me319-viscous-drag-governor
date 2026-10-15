# A Viscous-Drag Governor (VDG)

**ME319: Control Systems (Autumn 2026) — Assignment 1: Reinventing the Steam Engine Governor**

An automatic steam-engine speed regulator that senses shaft speed by **viscous shear in an
oil film** instead of by centrifugal force. A drag-cup rotor keyed to the engine shaft shears
an oil film that drags a concentric, torsion-spring-restrained sensing drum; the drum
deflection is therefore **linear in speed**, unlike the $\omega^2$ lift of a flyball governor.
The drum positions the steam throttle through a lever and link.

| | |
|---|---|
| Paper (PDF) | [`paper/main.pdf`](paper/main.pdf) — 6 pages + 1-page AI appendix |
| LaTeX source | [`paper/main.tex`](paper/main.tex) |
| Simulation | [`sim/vdg_main.m`](sim/vdg_main.m) (MATLAB / Octave) |
| Figures | [`figs/`](figs) |
| Numerical results | [`docs/results.txt`](docs/results.txt) |

---

## Governing equations

Deviation variables about the normalised operating point $\omega_0 = 1$, $T_{L0} = 1$:

```
J  d(δω)/dt                     = -b δω + Kt δxv - δTL
Js d²(δθ)/dt² + cs d(δθ)/dt + ks δθ = cs δω
tau_v d(δxv)/dt + δxv           = -Kl δθ
```

with transfer functions

```
Gp(s) = Kt / [ (J s + b)(tau_v s + 1) ]  =  1 / (s² + 10.1 s + 1)
Gc(s) = Kl cs / (Js s² + cs s + ks)      =  16 / (0.02 s² + 0.2 s + 1)
L(s)  = Gc(s) Gp(s)                      =  16 / (0.02 s⁴ + 0.402 s³ + 3.04 s² + 10.3 s + 1)
G_CL(s) = L / (1 + L)
```

(The sensor and linkage values are quoted after dividing the drum equation by `ks`, so the
numbers that set the dynamics are the ratios `cs/ks`, `Js/ks` and `Kl cs/ks`.)

Adding a viscous reset dashpot (time constant $T_i$) gives

```
Gc(s) = [Kl cs / (Js s² + cs s + ks)] · (Ti s + 1)/(Ti s)
```

## Parameters

| Symbol | Value | Meaning |
|---|---|---|
| `J` | 10 | engine + flywheel inertia (≈ 5.1×10³ kg·m² at 60 rpm / 20 kW) |
| `b` | 1 | effective viscous damping, `b·ω0 = T_L0 = 1` |
| `Kt` | 1 | throttle-to-torque gain — the nominal balance `Kt·xv0 = T_L0 = 1` |
| `tau_v` | 0.1 s | throttle-valve time constant |
| `Js` | 4×10⁻⁴ | sensing-drum inertia (≈ 0.20 kg·m²; `Js/ks = 0.02 s²`) |
| `cs` | 4×10⁻³ | oil-film shear coefficient (≈ 2.0 N·m·s/rad; transmits 12.7 N·m, 80 W) |
| `ks` | 0.02 | torsion-spring stiffness (θ_ss = c_s ω0/k_s: 72° at rated speed) |
| `Kl` | 80 | lever/linkage gain (0.4 m drum radius feeding a 20:1 compound lever) |
| `Ti` | 4 s | reset dashpot time constant (PI variant) |

Nominal balance: `Kt·xv0 = T_L0 = 1`. The damping torque at rated speed is *inside* `TL`, so
it must not be counted a second time in the steady-state balance; an earlier revision used
`Kt = 2`, which double-counted it. The sensor and linkage values were fitted (see the
appendix of the paper for exactly what was fitted and what was derived).

## Key results

* Closed-loop poles: `-2.432 ± 2.685j`, `-7.618 ± 2.595j` rad/s
* Dominant pair: `ωn = 3.62 rad/s`, `ζ = 0.671`
* Routh first column: `0.02, 0.402, 2.528, 7.596, 17` — all positive ⇒ asymptotically stable
* Phase margin 66.3°, gain margin 12.0 dB
* Steady-state droop `1/(1+G0) = 5.88 %` per unit load, i.e. 1.76 % for a 30 % load step

| Case (+30 % load step) | ω∞ | error | peak dev. | t_s | x_v∞ |
|---|---|---|---|---|---|
| Open loop | 0.700 | 30.0 % | 30.0 % | 38.96 s | 1.000 |
| VDG, proportional | 0.982 | 1.76 % | 1.85 % | 1.58 s (±2 %) | 1.282 |
| VDG, + reset dashpot | 1.000 | 0.00 % | 1.77 % | 5.5 s (±0.5 %) | 1.300 |

---

## Reproducing the results

**Requirements.** MATLAB R2019b or newer to regenerate the figures and numbers; a LaTeX
distribution with the `IEEEtran` class to rebuild the PDF. No toolbox functions are used —
the time response is integrated with an explicit fixed-step RK4 scheme and the poles, Routh
array and margins are computed from the polynomial coefficients directly, so the scripts
should also run under GNU Octave. **The scripts were validated in MATLAB R2025b**; Octave was
not available in the development environment, so Octave compatibility is by construction and
has not been tested.

**1. Regenerate every number and figure**

```bash
cd sim
matlab -batch "vdg_main"          # MATLAB
octave --eval "vdg_main"          # GNU Octave (alternative)
```

This writes

* `../figs/fig_open_closed.pdf|png`, `fig_reset`, `fig_poles`, `fig_zoom`, `fig_bode`
* `../docs/results.txt` — the full numerical log (transfer functions, poles, Routh first
  column, stability margins, transient metrics)

and prints the same log to the console.

**2. Rebuild the paper**

```bash
cd ../paper
pdflatex main.tex && pdflatex main.tex      # twice, for cross-references
```

The PDF is written to `paper/main.pdf`. The LaTeX source picks up the figures from `../figs/`.

### Source files

| File | Purpose |
|---|---|
| `sim/vdg_params.m` | parameter table with the physical justification of each value |
| `sim/vdg_model.m` | state-space + transfer-function model (plant, controller, loop, closed loop) |
| `sim/vdg_routh.m` | Routh array and stability verdict |
| `sim/vdg_rk4.m` | fixed-step RK4 integrator for the state-space |
| `sim/vdg_main.m` | driver: analysis, simulation, figures, results table |

## Limitations of the model

The nominal model assumes rigid, massless, backlash-free linkages, an isothermal oil film and
small perturbations about the operating point. Oil-viscosity drift with temperature, pivot
friction, lever backlash and spring creep are the dominant physical limitations; they are
analysed in Sec. VI-C of the paper.

## License

MIT — see [`LICENSE`](LICENSE).
