# SphereTurbIC simulation setup

This setup follows the collapse of a spherical cloud with an initial turbulent
velocity field. Self-gravity, hydrodynamics, single-group radiation transport,
gas-dust thermal coupling, tabulated cooling, and Sink particles are enabled
from the start. There is no time-dependent turbulent forcing.

- Problem generator: [testSphereTurbIC.cpp](../src/problems/SphereTurbIC/testSphereTurbIC.cpp)
- Runtime input: [SphereTurbIC_balanced_2tff.toml](SphereTurbIC_balanced_2tff.toml)

## Cloud and ambient medium

| Quantity | Initial value |
| --- | --- |
| Analytic cloud mass | 1000 solar masses |
| Cloud number density | 100 cm^-3 |
| Ambient number density | 1 cm^-3 |
| Cloud radius, R | 3.460368 pc |
| Cubic box side, L | 27.682946 pc = 8R = 4 cloud diameters |
| Cloud gas temperature | 48.9 K, approximately 49 K |
| Ambient gas temperature | Approximately 4900 K, determined by the EOS at the cloud pressure |

Here the reported number density is defined by

$$
\rho = 2.33\,m_{\rm H}\,n,
\qquad m_{\rm H}=1.6735575\times10^{-24}\ {\rm g}.
$$

Thus the cloud density is 3.899388975e-22 g cm^-3, and its radius is set by

$$
R=\left(\frac{3M}{4\pi\rho_{\rm cloud}}\right)^{1/3}.
$$

The PGEN computes the cloud pressure with the tabulated EOS and then sets the
ambient internal energy at that same pressure. This is **initial gas-pressure
balance**, not a claim of thermal-source equilibrium. Cooling, heating, and
radiation can subsequently change both temperatures.

## Initial turbulence: decaying, not continuously driven

Generate the velocity cube using the repository's perturbation script:

```sh
cd src/turbulence
python3 perturbation.py --kmin=2 --kmax=64 --size=128 --alpha=2 --f_solenoidal=0.5
```

The spectral parameter `--alpha=2` specifies an initial shell-integrated velocity
energy spectrum approximately proportional to k^-2, not a Kolmogorov k^-5/3
spectrum. The dimensionless wave number is k = L/lambda. With kmin = 2, the
longest input wavelength is L/2 = 13.841473 pc; the k = 1 box fundamental is not
included. `--f_solenoidal=0.5` selects a mixture of solenoidal and compressive
components. This spectral `alpha` is separate from the virial parameter below.

The PGEN subtracts the mean velocity measured inside the initial cloud and
normalizes the cloud's three-dimensional velocity dispersion to

$$
\alpha_{\rm kin}
=\frac{5\sigma_{\rm 3D}^{2}R}{3GM}=1,
\qquad
\sigma_{\rm 3D}=\left(\frac{3GM}{5R}\right)^{1/2}
\simeq0.864\ {\rm km\,s^{-1}}.
$$

## Mesh and Jeans refinement

The base grid is 128^3, with `amr.max_level = 7` and a factor of two between
successive levels. The finest **allowed** spacing is

$$
\Delta x_{\min} = \frac{L}{128\,2^7}
\simeq348.5\ {\rm AU}.
$$

This corresponds to an effective 16384^3 grid only where level 7 is actually
created; the whole domain does not start at that resolution.

The PGEN tags cells for refinement when both

$$
n\ge10\ {\rm cm^{-3}},
\qquad
\lambda_J=c_s\sqrt{\frac{\pi}{G\rho}}<10\,\Delta x
$$

hold, until the maximum allowed level is reached. Here c_s is the local
isothermal sound speed returned by the EOS. The AMR Jeans resolution is N = 10;
the Sink creation threshold uses N = 4 instead.

## Radiation, cooling, and gas-dust coupling

- The initial radiation energy density is spatially uniform,
  E_rad = a_rad (10 K)^4, with zero initial radiation flux.
- Radiation is evolved dynamically in one grey group. The initial radiation
  temperature is not a prescribed, permanently fixed dust temperature or an
  externally maintained ISRF.
- The PGEN sets constant mass opacities
  kappa_P = kappa_E = kappa_F = 2 cm^2 g^-1 of gas.
- The reduced speed of light is c_hat/c = 1e-4. The radiation-energy floor
  corresponds to a blackbody temperature of 2.725 K.
- `radiation.dust_gas_interaction_coeff = 2.5e-32`, in the code's units of
  erg cm^3 s^-1 K^-3/2. This is 100 times the default 2.5e-34. The separate
  dust-fluid dynamics module is disabled; gas-dust **thermal coupling** is enabled.
- The cooling table is
  `extern/cooling/CloudyData_UVB=HM2012_shielded_resampled_noPE.h5`.
  This is the shielded HM2012 table without photoelectric heating.
- The runtime gas temperature floor is 10 K.

## Boundary conditions

- **Gas:** reflecting on all six box faces.
- **Radiation:** `amrex::BCType::foextrap` for radiation energy and every
  radiation-flux component on all six faces. These per-component overrides
  are set in the PGEN; the gas retains its reflecting conditions.
- **Gravity:** isolated/open gravity through the OpenBC solver for the fully
  non-periodic domain; gas reflecting boundaries do not imply reflecting gravity.

## Sink creation and accretion floor

For this non-MHD problem the Jeans threshold is

$$
\rho>\rho_J=\frac{\pi c_s^2}{G(4\,\Delta x)^2},
\qquad\text{equivalently}\qquad\lambda_J<4\,\Delta x.
$$

The PGEN also specifies a fixed **remaining-gas accretion density floor**:

$$
\rho_{\rm acc,floor}
=2.33\,m_{\rm H}\,(100\ {\rm cm^{-3}})
=3.899388975\times10^{-22}\ {\rm g\,cm^{-3}}.
$$

## Duration and outputs

The initial cloud free-fall time is

$$
t_{\rm ff,0}=\sqrt{\frac{3\pi}{32G\rho_{\rm cloud}}}
=3.370980\ {\rm Myr}.
$$

The configured physical end time is 2 initial free-fall times, or 6.741959 Myr
(`stop_time = 212760047867775.34` seconds). This describes the requested run
duration, not a claim that a particular job has reached it. Plotfiles are
requested every 0.01 Myr and checkpoints every 0.1 Myr.
