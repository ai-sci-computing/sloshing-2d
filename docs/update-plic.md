## Update: Dambreak with Baffle

A dambreak comparison against OpenFOAM's interFoam revealed a weakness: the wavefront was advancing at 28% of Ritter's analytical speed. Published 2D CFD lands between 70% and 95%.

### Changes

**Velocity extrapolation.** The semi-Lagrangian back-trace routinely lands in the air region above the free surface, where nothing sets the velocity. Every back-trace into air returns zero, damping surface motion each step. Fix: propagate fluid velocities into the air band before each advection step, weighted by the VOF gradient along the interface normal (Foster & Fedkiw 2001).

**Geometric VOF (PLIC).** The algebraic VOF scheme (Superbee TVD) smeared air into the water bulk. During oscillation, surface cells lose VOF to numerical diffusion; on flow reversal, depleted values ratchet inward. This is structural to diffuse-interface tracking — no limiter tuning fixes it. CICSAM eliminated the smearing but introduced diagonal flotsam, a known NVD artefact. Fix: PLIC — reconstruct a straight line segment per interface cell from a Youngs-estimated normal, compute fluxes as exact polygon intersections. Bulk cells stay at exactly 0 or 1 by construction.

**Boundary conditions.** The SOLID border layer and explicit `enforce_boundary_conditions()` were replaced by folding wall BCs into the pressure Poisson equation as non-homogeneous Neumann. This exposed a CSF surface tension bug near internal obstacles (baffle VOF = 0 was being read as a fake water-air interface).

### Numbers

| metric                         | before | after | change   |
| ------------------------------ | ------ | ----- | -------- |
| Sloshing frequency error       | 5.9%   | 0.8%  | -86%     |
| Ritter wavefront ratio         | 0.28   | 0.71  | +154%    |
| Thin-film spreading ratio      | --     | 0.64  | new test |
| Hydrostatic RMS error          | 5.1%   | 5.1%  | =        |
| Volume drift (violent shaking) | ~0%    | 0%    | =        |

The Ritter ratio at 0.71 sits inside the published range for 2D CFD. Volume conservation is exact to machine precision with PLIC.
