#include "SpaceHub/src/spaceHub.hpp"
#include "SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
using namespace orbit;
using f = Interactions<NewtonianGrav>;
using Solver = methods::DefaultMethod<f>;
using Particle = Solver::Particle;

void job(double inc_in) {
    std::string file_name = "system-check-" + std::to_string(int(inc_in)) + ".txt";

    Particle HD3167{0.94_Ms};
    // Particle HD3167b{5.02_Me};
    Particle HD3167d{6.9_Me};
    Particle HD3167c{9.8_Me};

    double mutual_angle = inc_in * 1_deg;

    auto d_orb = Elliptic(HD3167.mass, HD3167d.mass, 0.07757, 0.0, 89.3_deg + mutual_angle, 0.0, 0.0, isotherm);

    auto c_orb = Elliptic(HD3167.mass, HD3167c.mass, 0.1841, 0.05, 89.3_deg, 0.0, 0.0, isotherm);

    move_particles(d_orb, HD3167d);
    move_particles(c_orb, HD3167c);

    move_to_COM_frame(HD3167, HD3167d, HD3167c);

    Solver sim{0, HD3167, HD3167d, HD3167c};

    Solver::RunArgs args;

    args.atol = 1e-13;

    double end_time_limit = 1e8_year;

    args.add_stop_condition(end_time_limit);

    /*int ejection_id = 0;

    auto any_ejection = [&](auto& ptc, auto h) {
        auto [a_d, e_d] = calc_a_e(ptc.mass(0) + ptc.mass(1), ptc.pos(0) - ptc.pos(1), ptc.vel(0) - ptc.vel(1));
        auto [a_c, e_c] = calc_a_e(ptc.mass(0) + ptc.mass(2), ptc.pos(0) - ptc.pos(2), ptc.vel(0) - ptc.vel(2));
        auto [a_x, e_x] = calc_a_e(ptc.mass(0) + ptc.mass(3), ptc.pos(0) - ptc.pos(3), ptc.vel(0) - ptc.vel(3));

        if (e_d >= 1.0) {
            ejection_id = 1;
            return true;
        } else if (e_c >= 1.0) {
            ejection_id = 2;
            return true;
        } else if (e_x >= 1.0) {
            ejection_id = 3;
            return true;
        } else {
            return false;
        }
    };
    args.add_stop_condition(StepSlice(any_ejection, 10));*/

    args.add_operation(StepSlice(DefaultWriter(file_name), 50));

    sim.run(args);
}

int main(int argc, char** argv) {
    tf::Executor executor;

    executor.silent_async(job, 0.0);

    executor.silent_async(job, 5.0);

    executor.silent_async(job, 10.0);

    executor.silent_async(job, 15.0);

    executor.silent_async(job, 20.0);

    executor.silent_async(job, 25.0);

    executor.silent_async(job, 30.0);

    executor.silent_async(job, 35.0);

    executor.silent_async(job, 40.0);

    executor.silent_async(job, 45.0);

    executor.silent_async(job, 50.0);

    executor.silent_async(job, 55.0);

    executor.silent_async(job, 60.0);

    executor.silent_async(job, 65.0);

    executor.silent_async(job, 70.0);

    executor.silent_async(job, 75.0);

    executor.silent_async(job, 80.0);

    executor.silent_async(job, 85.0);

    executor.silent_async(job, 90.0);

    executor.wait_for_all();

    return 0;
}
