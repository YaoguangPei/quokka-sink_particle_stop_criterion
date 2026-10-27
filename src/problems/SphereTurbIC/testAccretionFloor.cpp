#include "QuokkaSimulation.hpp"
#include "particles/particle_accretion.hpp"
#include <cmath>

struct FixedAccretionFloorTest {};
template <> struct quokka::EOS_Traits<FixedAccretionFloorTest> {
	static constexpr double gamma = 5. / 3.;
	static constexpr double mean_molecular_weight = C::m_u;
};
template <> struct Physics_Traits<FixedAccretionFloorTest> : DefaultPhysicsTraits {
	static constexpr bool is_hydro_enabled = true;
	static constexpr bool is_self_gravity_enabled = true;
};
template <> struct Particle_Traits<FixedAccretionFloorTest> : DefaultParticleTraits {
	static constexpr double sink_accretion_density_floor = 1.;
	static constexpr ParticleSwitch particle_switch = ParticleSwitch::Sink;
};

auto problem_main() -> int
{
	init_extern_parameters();
	amrex::Real small_temperature = 1.e-10;
	amrex::Real small_density = 1.e-100;
	eos_init(small_temperature, small_density);
	const amrex::Box box(amrex::IntVect(0), amrex::IntVect(3));
	const amrex::BoxArray boxes(box);
	const amrex::DistributionMapping mapping(boxes);
	const amrex::RealBox physical({0., 0., 0.}, {4., 4., 4.});
	const int periodic[3] = {0, 0, 0};
	const amrex::Geometry geometry(box, &physical, 0, periodic);
	amrex::MultiFab state(boxes, mapping, Physics_Indices<FixedAccretionFloorTest>::nvarTotal_cc, 0);
	amrex::MultiFab request(boxes, mapping, 1, 0), scale(boxes, mapping, 1, 0), floor(boxes, mapping, 1, 0);
	floor.setVal(0.001);
	for (double rho : {0.5, 1., 1.1, 2.}) {
		state.setVal(0.);
		state.setVal(rho, HydroSystem<FixedAccretionFloorTest>::density_index, 1, 0);
		state.setVal(1.e12 * rho / (2. / 3.), HydroSystem<FixedAccretionFloorTest>::internalEnergy_index, 1, 0);
		state.setVal(1.e12 * rho / (2. / 3.), HydroSystem<FixedAccretionFloorTest>::energy_index, 1, 0);
		request.setVal(-0.9);
		scale.setVal(1.);
		quokka::SinkAccretionUtils::ComputeScaleDown<FixedAccretionFloorTest>(state, request, scale, geometry, nullptr, floor);
		const double expected = rho <= 1. ? 0. : std::max(-0.25, 1. / rho - 1.);
		const double low = request.min(0), high = request.max(0);
		AMREX_ALWAYS_ASSERT(std::abs(low - expected) < 1.e-12 && std::abs(high - expected) < 1.e-12);
		AMREX_ALWAYS_ASSERT((1. + low) * rho >= std::min(rho, 1.) - 1.e-12);
		amrex::Print() << "ACCRETION_FLOOR_TEST rho=" << rho << " remaining=" << (1. + low) * rho << " passed\n";
	}
	return 0;
}
