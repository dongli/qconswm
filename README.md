Introduction
============

This is an shallow-water (i.e. barotropic) implemention of the
quadratic-conservation dynamical core (dycore for short) developed in LASG lab
at institute of atmospheric physics (IAP) by Bin Wang's team. This type of
dycore can preserve total mass, total energy and some other conservation
quantities exactly to machine precision level.

The total energy or quadratic form of IAP transformed variables is linked with
computational stability, so when total energy is conserved, the long-term
computing can be ensured which is key to climate simulation.

Currently, the implemention contains Arakawa C grid, second-order spatial
difference satisfying antisymmetric condition, and second-order
predict-correct time integrator, etc. More options will be added soon.

| Time integrator | Order | Implicit | Available |
|-----------------|-------|----------|-----------|
| predict-correct | 2nd   | No       | Yes       |
| middle-point    | 2nd   | Yes      | Yes       |
| runge-kutta     | xxx   | Yes      | Not yet   |
| leap-frog       | 3rd   | Yes      | Not yet   |

Authors
=======

- Bin Wang
- Li Dong
- ...
