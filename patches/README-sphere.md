# SphereTurbIC review snapshot

This snapshot uses upstream Quokka commit `1e792e10295a22e8876e29cbb33edec584c775e1`
plus the radiation changes actually integrated on SilkRiver from PR2326 and
PR2344. PR2355 (the continuous turbulence amplitude controller) is excluded.
The velocity field is an initial condition; `turbulence.enabled = 0`.

Main problem: `src/problems/SphereTurbIC/`.
Main input: `input/SphereTurbIC_balanced_2tff.toml`, with stop_time = 2 tff.
The original misplaced root-level PGEN and TOML are not present in this tree.

The AMReX dependency remains pinned to the tested base commit
`1de18774af5d8a368fa2a398f048e205a2c9f84f`. Before configuring, initialize the
submodules and apply the exact locally integrated OpenBC fix derived from
AMReX PR5664 (not a claim that every file equals the latest PR head):

```sh
git submodule update --init --recursive
git -C extern/amrex apply --check ../../patches/amrex-openbc-5664.patch
git -C extern/amrex apply ../../patches/amrex-openbc-5664.patch
```

The current PGEN's optional uniform volumetric heating requires the accompanying
QuokkaSimulation.hpp and ResampledCooling.hpp changes included here. The fixed
remaining-gas Sink accretion floor and its supporting particle changes are also
included. This is the later thermally balanced setup; the figures from job167827
were produced with the older initial conditions. This is not a noPE run: the
original cooling table, RT, and 100x gas-dust thermal coupling are retained.
The ambient equilibrium branch is thermally unstable to isobaric perturbations.

Run from the repository root, with the same input data at:
`extern/cooling/CloudyData_UVB=HM2012_shielded_resampled.h5` and
`src/turbulence/zdrv.hdf5`. The input uses relative paths for review portability.
The requested output intervals remain 0.01 Myr (plots) and 0.1 Myr (checkpoints).
No production job, source working tree, or server input was modified by this export.
