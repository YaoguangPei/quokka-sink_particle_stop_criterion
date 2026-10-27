#define problem_main sphere_production_problem_main
#include "testSphereTurbIC.cpp"
#undef problem_main
#include <fstream>

namespace {
void writePhysicalVelocityReference()
{
    std::string filename;
    amrex::ParmParse("perturb").get("filename", filename);
    const hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    AMREX_ALWAYS_ASSERT(file >= 0);
    constexpr int N = 128;
    std::array<std::vector<double>, 3> raw;
    const char *names[3] = {"pertx", "perty", "pertz"};
    for (int d = 0; d < 3; ++d) {
        const hid_t dataset = H5Dopen2(file, names[d], H5P_DEFAULT);
        const hid_t space = H5Dget_space(dataset);
        hsize_t dims[3]{};
        AMREX_ALWAYS_ASSERT(H5Sget_simple_extent_dims(space, dims, nullptr) == 3);
        AMREX_ALWAYS_ASSERT(dims[0] == N && dims[1] == N && dims[2] == N);
        raw[d].resize(N*N*N);
        AMREX_ALWAYS_ASSERT(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw[d].data()) >= 0);
        H5Sclose(space); H5Dclose(dataset);
    }
    H5Fclose(file);
    std::vector<amrex::GpuArray<int, 3>> cells;
    amrex::GpuArray<double, 3> mean{};
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                const double r2 = std::pow(i+.5-64.,2)+std::pow(j+.5-64.,2)+std::pow(k+.5-64.,2);
                if (r2 > 64.) { continue; }
                cells.push_back({i,j,k});
                for (int d=0;d<3;++d) { mean[d] += raw[d][(i*N+j)*N+k]; }
            }
        }
    }
    for (double &value : mean) { value /= static_cast<double>(cells.size()); }
    double variance=0.;
    for (auto const &q : cells) {
        for (int d=0;d<3;++d) { variance += std::pow(raw[d][(q[0]*N+q[1])*N+q[2]]-mean[d],2); }
    }
    variance /= static_cast<double>(cells.size());
    const double rho = cloud_number_density * particle_mass;
    const double rad = std::cbrt(3.*cloud_mass_msun*C::M_solar/(4.*M_PI*rho));
    const double target = std::sqrt((3./5.)*C::Gconst*cloud_mass_msun*C::M_solar/rad);
    std::ofstream out("velocity_reference.csv"); out << std::setprecision(17) << "i,j,k,vx_kms,vy_kms,vz_kms\n";
    for (auto const &q : cells) {
        out << q[0] << ',' << q[1] << ',' << q[2];
        for (int d=0;d<3;++d) { out << ',' << (raw[d][(q[0]*N+q[1])*N+q[2]]-mean[d])*target/std::sqrt(variance)/1.e5; }
        out << '\n';
    }
    out.close(); AMREX_ALWAYS_ASSERT(out.good());
    amrex::Print() << "RAW_HDF5_VELOCITY_REFERENCE_WRITTEN cells=" << cells.size() << '\n';
}

void checkHdf5Layout()
{
    constexpr int N = 8;
    const char *filename = "velocity_layout_unit.h5";
    const hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    AMREX_ALWAYS_ASSERT(file >= 0);
    const hsize_t dimensions[3] = {N, N, N};
    const hid_t space = H5Screate_simple(3, dimensions, nullptr);
    const double coefficients[3][3] = {{3., 11., 37.}, {7., -13., 2.}, {-5., 17., 29.}};
    const char *names[3] = {"pertx", "perty", "pertz"};
    for (int d = 0; d < 3; ++d) {
        std::vector<double> values(N * N * N);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    values[(i * N + j) * N + k] = coefficients[d][0] * i + coefficients[d][1] * j + coefficients[d][2] * k;
                }
            }
        }
        const hid_t dataset = H5Dcreate2(file, names[d], H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        AMREX_ALWAYS_ASSERT(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()) >= 0);
        H5Dclose(dataset);
    }
    H5Sclose(space); H5Fclose(file);
    turb_data data;
    std::string input(filename);
    initialize_turbdata(data, input);
    auto tx = get_tabledata(data.dvx); auto ty = get_tabledata(data.dvy); auto tz = get_tabledata(data.dvz);
    const amrex::GpuArray<int, 3> size{N, N, N};
    const amrex::GpuArray<double, 3> lo{0., 0., 0.};
    const amrex::GpuArray<double, 3> length{N, N, N};
    for (auto const &position : {amrex::GpuArray<double, 3>{1.5, 2.5, 4.5}, amrex::GpuArray<double, 3>{2.25, 3.1, 1.8}}) {
        const double measured[3] = {sample(tx.const_table(), size, position, lo, length),
                                   sample(ty.const_table(), size, position, lo, length),
                                   sample(tz.const_table(), size, position, lo, length)};
        for (int d = 0; d < 3; ++d) {
            double expected = 0.;
            for (int axis = 0; axis < 3; ++axis) { expected += coefficients[d][axis] * (position[axis] - .5); }
            AMREX_ALWAYS_ASSERT(std::abs(measured[d] - expected) < 1.e-11);
        }
    }
    amrex::Print() << "HDF5_LAYOUT_UNIT_PASSED integer_and_fractional_physical_coordinates\n";
}

void coolingStep(double rho, double &energy, double dt, QuokkaSimulation<SphereTurbIC> const &sim)
{
    quokka::valarray<double, 1> y{energy};
    const double floor = quokka::EOS<SphereTurbIC>::ComputeEintFromTgas(rho, 10.);
    const quokka::valarray<double, 1> atol{.01 * floor};
    int substeps = 0;
    const quokka::ResampledCooling::ResampledCoolingFunctor rhs(rho, sim.resampledTables_.const_tables_host(), background_heating_rate);
    rk_adaptive_integrate(rhs, 0., y, dt, 1.e-4, atol, substeps);
    AMREX_ALWAYS_ASSERT(substeps < maxStepsODEIntegrate && std::isfinite(y[0]) && y[0] > 0.);
    energy = y[0];
}

void radiationStep(double rho, double &energy, double &erad, double dt, QuokkaSimulation<SphereTurbIC> const &sim)
{
    CouplingCell<SphereTurbIC> cell{};
    cell.Egas0 = energy; cell.Erad0[0] = erad; cell.rho = rho; cell.dt = dt;
    cell.tau_scale = dt * RadSystem<SphereTurbIC>::c_hat_;
    const double nH = RadSystem<SphereTurbIC>::ComputeNumberDensityH(rho, cell.massScalars);
    cell.dtK = dt * sim.dustGasInteractionCoeff_ * nH * nH;
    cell.Tfloor = 10.; cell.Emin = quokka::EOS<SphereTurbIC>::ComputeEintFromTgas(rho, 10.);
    const auto result = RadSystem<SphereTurbIC>::SolveDustCoupling(cell, 1.e-12);
    AMREX_ALWAYS_ASSERT(result.converged && std::isfinite(result.Egas));
    const double before = energy + 1.e4 * erad;
    energy = result.Egas; erad = result.Erad[0];
    AMREX_ALWAYS_ASSERT(std::abs((energy + 1.e4 * erad) / before - 1.) < 1.e-12);
}
} // namespace

auto problem_main() -> int
{
    AMREX_ALWAYS_ASSERT(amrex::ParallelDescriptor::NProcs() == 1);
    checkHdf5Layout();
    writePhysicalVelocityReference();
    balance_thermal_sources = 1;
    auto bc = quokka::BC<SphereTurbIC>(quokka::BCType::reflecting);
    QuokkaSimulation<SphereTurbIC> sim(bc);
    establishThermalBalance(sim);
    const double rho_c = cloud_number_density * particle_mass;
    const double initial_e = quokka::EOS<SphereTurbIC>::ComputeEintFromTgas(rho_c, cloud_temperature);
    const double pressure = quokka::EOS<SphereTurbIC>::ComputePressure(rho_c, initial_e);
    const double tff = std::sqrt(3. * M_PI / (32. * C::Gconst * rho_c));
    for (double n : {cloud_number_density, ambient_number_density}) {
        const double rho = n * particle_mass;
        AMREX_ALWAYS_ASSERT(std::abs(initialGasRate(rho, pressure, sim) + background_heating_rate) < 1.e-4 * background_heating_rate);
        // Exercise the actual simulation -> cooling-integrator path, not only its scalar functor.
        const amrex::Box box(amrex::IntVect(0), amrex::IntVect(3));
        const amrex::BoxArray boxes(box);
        const amrex::DistributionMapping mapping(boxes);
        amrex::MultiFab state(boxes, mapping, Physics_Indices<SphereTurbIC>::nvarTotal_cc, 0);
        std::array<amrex::MultiFab, AMREX_SPACEDIM> empty_fc;
        const double test_e = quokka::EOS<SphereTurbIC>::ComputeEintFromPres(rho, pressure);
        state.setVal(0.);
        state.setVal(rho, HydroSystem<SphereTurbIC>::density_index, 1, 0);
        state.setVal(test_e, HydroSystem<SphereTurbIC>::internalEnergy_index, 1, 0);
        state.setVal(test_e, HydroSystem<SphereTurbIC>::energy_index, 1, 0);
        const double short_dt = 1.e6;
        const bool success = sim.addStrangSplitSourcesWithBuiltin<QuokkaSimulation<SphereTurbIC>::SourceOrder::forward>(
            state, empty_fc, nullptr, 0, 0., short_dt);
        const double observed = (state.min(HydroSystem<SphereTurbIC>::internalEnergy_index) - test_e) / short_dt;
        const double expected = quokka::ResampledCooling::resampled_cooling_function(rho, test_e, sim.resampledTables_.const_tables_host()) + background_heating_rate;
        AMREX_ALWAYS_ASSERT(success && std::abs(observed - expected) < 1.e-3 * background_heating_rate);
        amrex::Print() << "VOLUMETRIC_HEATING_PIPELINE_PASSED n=" << n << " measured_rate=" << observed << " expected_rate=" << expected << '\n';
        const double rho_lo = .999 * rho, rho_hi = 1.001 * rho;
        const double e_lo = quokka::EOS<SphereTurbIC>::ComputeEintFromPres(rho_lo, pressure);
        const double e_hi = quokka::EOS<SphereTurbIC>::ComputeEintFromPres(rho_hi, pressure);
        const double T_lo = quokka::EOS<SphereTurbIC>::ComputeTgasFromEint(rho_lo, e_lo);
        const double T_hi = quokka::EOS<SphereTurbIC>::ComputeTgasFromEint(rho_hi, e_hi);
        const double slope = (initialGasRate(rho_lo, pressure, sim) - initialGasRate(rho_hi, pressure, sim)) / (T_lo - T_hi);
        const double cp = (5. / 3.) * (e_lo / rho_lo - e_hi / rho_hi) / (T_lo - T_hi);
        amrex::Print() << "ISOBARIC_THERMAL_RESPONSE n=" << n << " dQdT=" << slope
            << " stable=" << (slope < 0.) << " linear_response_time_tff=" << rho * cp / std::abs(slope) / tff << '\n';
        for (double step_myr : {.01, .005}) {
            double energy = quokka::EOS<SphereTurbIC>::ComputeEintFromPres(rho, pressure);
            const double erad0 = C::a_rad * std::pow(radiation_temperature, 4);
            double erad = erad0, time = 0., max_drift = 0.;
            while (time < .1 * tff) {
                const double dt = std::min(step_myr * seconds_per_myr, .1 * tff - time);
                coolingStep(rho, energy, .5 * dt, sim);
                radiationStep(rho, energy, erad, dt, sim);
                coolingStep(rho, energy, .5 * dt, sim);
                const double now = quokka::EOS<SphereTurbIC>::ComputePressure(rho, energy);
                max_drift = std::max(max_drift, std::abs(now / pressure - 1.));
                time += dt;
            }
            amrex::Print() << std::setprecision(17) << "THERMAL_SOURCE_UNIT n=" << n << " dt_Myr=" << step_myr
                << " duration_tff=0.1 max_pressure_drift=" << max_drift << " radiation_relative_drift=" << erad / erad0 - 1.
                << " final_T=" << quokka::EOS<SphereTurbIC>::ComputeTgasFromEint(rho, energy) << '\n';
            AMREX_ALWAYS_ASSERT(max_drift < .02);
        }
    }
    amrex::Print() << "THERMAL_BALANCE_UNIT_PASSED gas_pressure_stability_and_dynamic_radiation\n";
    return 0;
}
