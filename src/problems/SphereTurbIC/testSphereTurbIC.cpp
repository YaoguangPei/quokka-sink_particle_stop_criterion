/// Self-gravitating sphere with a StarCluster-style turbulent velocity IC.
#include "QuokkaSimulation.hpp"
#include "fundamental_constants.H"
#include "hydro/hydro_system.hpp"
#include "radiation/radiation_system.hpp"
#include "turbulence/TurbDataReader.hpp"
#include "util/BC.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>

struct SphereTurbIC {};

namespace
{
constexpr double hydrogen_mass = 1.6735575e-24; // hydrogen atom mass [g], not atomic mass unit
constexpr double particle_mass = 2.33 * hydrogen_mass;
constexpr double seconds_per_myr = 3.15576e13;
double cloud_mass_msun = 1000.;
double cloud_number_density = 100.;
double ambient_number_density = 1.;
double cloud_temperature = 47.962492751223152;
double radiation_temperature = 10.;
int jeans_refinement_ncells = 10;
double radius = 0.;

template <typename Table>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE auto sample(Table const &field, amrex::GpuArray<int, 3> const &size, amrex::GpuArray<double, 3> const &position,
						     amrex::GpuArray<double, 3> const &lo, amrex::GpuArray<double, 3> const &length) -> double
{
	amrex::GpuArray<int, 3> index{};
	amrex::GpuArray<double, 3> fraction{};
	for (int d = 0; d < 3; ++d) {
		const double q = (position[d] - lo[d]) / length[d] * size[d] - 0.5;
		index[d] = static_cast<int>(std::floor(q));
		fraction[d] = q - index[d];
	}
	double value = 0.;
	for (int a = 0; a < 2; ++a) {
		for (int b = 0; b < 2; ++b) {
			for (int c = 0; c < 2; ++c) {
				const int i = (index[0] + a + size[0]) % size[0];
				const int j = (index[1] + b + size[1]) % size[1];
				const int k = (index[2] + c + size[2]) % size[2];
				// The official reader wraps a C-order HDF5 buffer in a Fortran-indexed Table3D.
				// Undo its spatial x/z reversal while leaving the vector component unchanged.
				value += field(k, j, i) * (a ? fraction[0] : 1. - fraction[0]) * (b ? fraction[1] : 1. - fraction[1]) *
					 (c ? fraction[2] : 1. - fraction[2]);
			}
		}
	}
	return value;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE auto dustOpacity(double rho, double temperature, double coefficient) -> double
{
	const double n = rho / particle_mass;
	const double x = amrex::min(amrex::max((n - 50.) / 40., 0.), 1.);
	const double t = amrex::min(amrex::max(temperature, 1.), 150.) / 10.;
	return coefficient * x * x * (3. - 2. * x) * t * t;
}
} // namespace

template <> struct quokka::EOS_Traits<SphereTurbIC> {
	static constexpr double gamma = 5. / 3.;
	static constexpr double mean_molecular_weight = C::m_u;
	using EOSBackend = quokka::EOSTabulated<SphereTurbIC>;
};
template <> struct HydroSystem_Traits<SphereTurbIC> {
	static constexpr bool reconstruct_eint = true;
};
template <> struct Physics_Traits<SphereTurbIC> : DefaultPhysicsTraits {
	static constexpr bool is_hydro_enabled = true;
	static constexpr bool is_self_gravity_enabled = true;
	static constexpr bool is_radiation_enabled = true;
	static constexpr bool is_dust_enabled = false;
};
template <> struct Particle_Traits<SphereTurbIC> : DefaultParticleTraits {
	static constexpr ParticleSwitch particle_switch = ParticleSwitch::Sink;
	static constexpr double sink_accretion_density_floor = 100. * particle_mass;
};
template <> struct RadSystem_Traits<SphereTurbIC> {
	static constexpr double c_hat_over_c = 1.e-4;
	static constexpr double Erad_floor = C::a_rad * 2.725 * 2.725 * 2.725 * 2.725;
	static constexpr int beta_order = 1;
};
template <> struct ISM_Traits<SphereTurbIC> {
	static constexpr bool enable_dust_gas_thermal_coupling_model = true;
};
template <> struct SimulationData<SphereTurbIC> {
	amrex::TableData<double, 3> dvx, dvy, dvz;
	amrex::GpuArray<int, 3> size{};
	amrex::GpuArray<double, 3> velocity_mean{};
	double velocity_scale = 0.;
};
template <> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE auto RadSystem<SphereTurbIC>::ComputePlanckOpacity(double rho, double temperature) -> amrex::Real
{
	return dustOpacity(rho, temperature, 0.1);
}
template <> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE auto RadSystem<SphereTurbIC>::ComputeFluxMeanOpacity(double rho, double temperature) -> amrex::Real
{
	return dustOpacity(rho, temperature, 0.0316);
}

template <> void QuokkaSimulation<SphereTurbIC>::preCalculateInitialConditions()
{
	amrex::ParmParse pp("perturb");
	std::string filename;
	double alpha = 1.;
	pp.get("filename", filename);
	pp.query("virial_parameter", alpha);
	AMREX_ALWAYS_ASSERT(alpha > 0.);
	turb_data data;
	initialize_turbdata(data, filename);
	auto vx = get_tabledata(data.dvx);
	auto vy = get_tabledata(data.dvy);
	auto vz = get_tabledata(data.dvz);
	auto const fx = vx.const_table();
	auto const fy = vy.const_table();
	auto const fz = vz.const_table();
	const auto low = vx.lo();
	const auto high = vx.hi();
	for (int d = 0; d < 3; ++d) {
		AMREX_ALWAYS_ASSERT(low[d] == 0);
		userData_.size[d] = high[d] + 1;
	}
	AMREX_ALWAYS_ASSERT(userData_.size[0] == userData_.size[1] && userData_.size[1] == userData_.size[2]);
	const auto lo = geom[0].ProbLoArray();
	const auto hi = geom[0].ProbHiArray();
	const auto dx = geom[0].CellSizeArray();
	const amrex::GpuArray<double, 3> length{hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
	const auto domain = geom[0].Domain();
	double count = 0., sum2 = 0.;
	amrex::GpuArray<double, 3> sum{};
	for (int k = domain.smallEnd(2); k <= domain.bigEnd(2); ++k) {
		for (int j = domain.smallEnd(1); j <= domain.bigEnd(1); ++j) {
			for (int i = domain.smallEnd(0); i <= domain.bigEnd(0); ++i) {
				const amrex::GpuArray<double, 3> pos{lo[0] + (i + .5) * dx[0], lo[1] + (j + .5) * dx[1], lo[2] + (k + .5) * dx[2]};
				double r2 = 0.;
				for (int d = 0; d < 3; ++d) {
					r2 += std::pow(pos[d] - .5 * (lo[d] + hi[d]), 2);
				}
				if (r2 > radius * radius) {
					continue;
				}
				const amrex::GpuArray<double, 3> v{sample(fx, userData_.size, pos, lo, length), sample(fy, userData_.size, pos, lo, length),
								   sample(fz, userData_.size, pos, lo, length)};
				for (int d = 0; d < 3; ++d) {
					sum[d] += v[d];
					sum2 += v[d] * v[d];
				}
				count += 1.;
			}
		}
	}
	AMREX_ALWAYS_ASSERT(count > 0.);
	double variance = sum2 / count;
	for (int d = 0; d < 3; ++d) {
		userData_.velocity_mean[d] = sum[d] / count;
		variance -= userData_.velocity_mean[d] * userData_.velocity_mean[d];
	}
	AMREX_ALWAYS_ASSERT(variance > 0. && std::isfinite(variance));
	const double sigma = std::sqrt(alpha * (3. / 5.) * Gconst_ * cloud_mass_msun * C::M_solar / radius);
	userData_.velocity_scale = sigma / std::sqrt(variance);
	userData_.dvx.resize(low, high);
	userData_.dvx.copy(vx);
	userData_.dvy.resize(low, high);
	userData_.dvy.copy(vy);
	userData_.dvz.resize(low, high);
	userData_.dvz.copy(vz);
	amrex::Print() << std::setprecision(17) << "TURB_IC alpha_kin_uniform_sphere=" << alpha << " sigma3D_kms=" << sigma / 1.e5
		       << " source_grid=" << userData_.size[0] << 'x' << userData_.size[1] << 'x' << userData_.size[2]
		       << " forcing=OFF normalization=cloud_base_cells_COM_subtracted\n";
}

template <> void QuokkaSimulation<SphereTurbIC>::setInitialConditionsOnGrid(quokka::grid const &grid)
{
	const auto state = grid.array_;
	const auto dx = grid.dx_;
	const auto lo = grid.prob_lo_;
	const auto hi = grid.prob_hi_;
	const amrex::GpuArray<double, 3> length{hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]};
	const double rho_cloud = cloud_number_density * particle_mass;
	const double rho_ambient = ambient_number_density * particle_mass;
	const double cloud_eint = quokka::EOS<SphereTurbIC>::ComputeEintFromTgas(rho_cloud, cloud_temperature);
	const double pressure = quokka::EOS<SphereTurbIC>::ComputePressure(rho_cloud, cloud_eint);
	const double ambient_eint = quokka::EOS<SphereTurbIC>::ComputeEintFromPres(rho_ambient, pressure);
	const double erad = C::a_rad * std::pow(radiation_temperature, 4);
	const double rcloud = radius;
	const auto fx = userData_.dvx.const_table();
	const auto fy = userData_.dvy.const_table();
	const auto fz = userData_.dvz.const_table();
	const auto size = userData_.size;
	const auto mean = userData_.velocity_mean;
	const double scale = userData_.velocity_scale;
	amrex::ParallelFor(grid.indexRange_, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
		const amrex::GpuArray<double, 3> pos{lo[0] + (i + .5) * dx[0], lo[1] + (j + .5) * dx[1], lo[2] + (k + .5) * dx[2]};
		double r2 = 0.;
		for (int d = 0; d < 3; ++d) {
			r2 += std::pow(pos[d] - .5 * (lo[d] + hi[d]), 2);
		}
		const bool inside = r2 <= rcloud * rcloud;
		const double rho = inside ? rho_cloud : rho_ambient;
		const double eint = inside ? cloud_eint : ambient_eint;
		const double weight = scale * rho / rho_cloud;
		const double vx = weight * (sample(fx, size, pos, lo, length) - mean[0]);
		const double vy = weight * (sample(fy, size, pos, lo, length) - mean[1]);
		const double vz = weight * (sample(fz, size, pos, lo, length) - mean[2]);
		for (int n = 0; n < Physics_Indices<SphereTurbIC>::nvarTotal_cc; ++n) {
			state(i, j, k, n) = 0.;
		}
		state(i, j, k, HydroSystem<SphereTurbIC>::density_index) = rho;
		state(i, j, k, HydroSystem<SphereTurbIC>::x1Momentum_index) = rho * vx;
		state(i, j, k, HydroSystem<SphereTurbIC>::x2Momentum_index) = rho * vy;
		state(i, j, k, HydroSystem<SphereTurbIC>::x3Momentum_index) = rho * vz;
		state(i, j, k, HydroSystem<SphereTurbIC>::internalEnergy_index) = eint;
		state(i, j, k, HydroSystem<SphereTurbIC>::energy_index) = eint + .5 * rho * (vx * vx + vy * vy + vz * vz);
		state(i, j, k, RadSystem<SphereTurbIC>::radEnergy_index) = erad;
	});
}

template <> void QuokkaSimulation<SphereTurbIC>::refineGrid(int lev, amrex::TagBoxArray &tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const double dx = geom[lev].CellSizeArray()[0];
	const double G = Gconst_;
	const int N = jeans_refinement_ncells;
	const auto state = state_new_cc_[lev].const_arrays();
	const auto tag = tags.arrays();
	amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k) noexcept {
		const double rho = state[bx](i, j, k, HydroSystem<SphereTurbIC>::density_index);
		const double cs = HydroSystem<SphereTurbIC>::ComputeIsothermalSoundSpeed(state[bx], i, j, k, nullptr);
		if (rho >= 10. * particle_mass && cs * std::sqrt(M_PI / (G * rho)) < N * dx) {
			tag[bx](i, j, k) = amrex::TagBox::SET;
		}
	});
}

template <>
void QuokkaSimulation<SphereTurbIC>::ComputeDerivedVar(int lev, std::string const &name, amrex::MultiFab &mf, int ncomp, amrex::MultiFab const &state_cc,
						       amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> const & /*state_fc*/) const
{
	const auto output = mf.arrays();
	if (name == "gpot") {
		const auto potential = phi[lev].const_arrays();
		amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE(int b, int i, int j, int k) noexcept { output[b](i, j, k, ncomp) = potential[b](i, j, k); });
		return;
	}
	AMREX_ALWAYS_ASSERT(name == "temperature" || name == "number_density" || name == "sound_speed" || name == "jeans_length");
	const int field = name == "temperature" ? 0 : (name == "number_density" ? 1 : (name == "sound_speed" ? 2 : 3));
	const double G = Gconst_;
	const auto state = state_cc.const_arrays();
	amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE(int b, int i, int j, int k) noexcept {
		const double rho = state[b](i, j, k, HydroSystem<SphereTurbIC>::density_index);
		const double eint = state[b](i, j, k, HydroSystem<SphereTurbIC>::internalEnergy_index);
		const double cs = HydroSystem<SphereTurbIC>::ComputeIsothermalSoundSpeed(state[b], i, j, k, nullptr);
		output[b](i, j, k, ncomp) = field == 0 ? quokka::EOS<SphereTurbIC>::ComputeTgasFromEint(rho, eint)
						       : (field == 1 ? rho / particle_mass : (field == 2 ? cs : cs * std::sqrt(M_PI / (G * rho))));
	});
}

auto problem_main() -> int
{
	amrex::ParmParse pp("problem");
	pp.query("cloud_mass_msun", cloud_mass_msun);
	pp.query("cloud_number_density", cloud_number_density);
	pp.query("ambient_number_density", ambient_number_density);
	pp.query("cloud_temperature", cloud_temperature);
	pp.query("radiation_temperature", radiation_temperature);
	pp.query("jeans_refinement_ncells", jeans_refinement_ncells);
	AMREX_ALWAYS_ASSERT(cloud_mass_msun > 0. && cloud_number_density > ambient_number_density && ambient_number_density > 0.);
	AMREX_ALWAYS_ASSERT(jeans_refinement_ncells >= 4);
	const double rho = cloud_number_density * particle_mass;
	radius = std::cbrt(3. * cloud_mass_msun * C::M_solar / (4. * M_PI * rho));
	const double tff = std::sqrt(3. * M_PI / (32. * C::Gconst * rho));
	auto bc = quokka::BC<SphereTurbIC>(quokka::BCType::reflecting);
	QuokkaSimulation<SphereTurbIC> sim(bc);
	const auto lo = sim.geom[0].ProbLoArray();
	const auto hi = sim.geom[0].ProbHiArray();
	for (int d = 0; d < 3; ++d) {
		AMREX_ALWAYS_ASSERT(std::abs((hi[d] - lo[d]) / (16. * radius) - 1.) < 1.e-10);
	}
	const double eint = quokka::EOS<SphereTurbIC>::ComputeEintFromTgas(rho, cloud_temperature);
	const double P = quokka::EOS<SphereTurbIC>::ComputePressure(rho, eint);
	const double rho_ambient = ambient_number_density * particle_mass;
	const double ambient_eint = quokka::EOS<SphereTurbIC>::ComputeEintFromPres(rho_ambient, P);
	amrex::Print() << std::setprecision(17) << "SPHERE_IC R_pc=" << radius / C::parsec << " box_pc=" << 16. * radius / C::parsec << " rho_cloud=" << rho
		       << " m_H_g=" << hydrogen_mass << " tff_Myr=" << tff / seconds_per_myr << " T_requested=" << cloud_temperature
		       << " T_EOS=" << quokka::EOS<SphereTurbIC>::ComputeTgasFromEint(rho, eint)
		       << " T_ambient=" << quokka::EOS<SphereTurbIC>::ComputeTgasFromEint(rho_ambient, ambient_eint) << " P_cloud=" << P
		       << " P_ambient=" << quokka::EOS<SphereTurbIC>::ComputePressure(rho_ambient, ambient_eint) << " refinement_N=" << jeans_refinement_ncells
		       << " sink_N=4 accretion_density_floor=" << Particle_Traits<SphereTurbIC>::sink_accretion_density_floor
		       << " self_gravity=ON RT_groups=1 forcing=OFF dust_gas_coeff=" << sim.dustGasInteractionCoeff_
		       << " automatic_thermal_balance=OFF extra_heating_vol=" << sim.volumetricHeatingRate_ << " velocity_HDF5_layout=C_order_corrected\n";
	sim.setInitialConditions();
	sim.evolve();
	return 0;
}
