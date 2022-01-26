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

void job(size_t thread_id, size_t scattering_num) {
    std::fstream file("stability-" + std::to_string(thread_id) + ".txt", std::ios::out);

    print(file, std::setprecision(16));

    print(std::cout, "starting job on thread ", thread_id, " \n");

    for (size_t i = 0; i < scattering_num; ++i) {
        double mass_x = pow(10.0, random::Uniform(-1, 1)) * 1_Mj;
        double a_x = 2 * pow(10.0, random::Uniform(-1, 1)) * 1_AU;

        Particle HD3167{0.94_Ms};
        // Particle HD3167b{5.02_Me};
        Particle HD3167d{6.9_Me};
        Particle HD3167c{9.8_Me};
        Particle HD3167x{mass_x};

        // auto p_orb1 = Elliptic(BH1.mass, BH2.mass, 0.0186, 0.0, consts::pi * double(retro), 0.0, 0.0, phi);

        double mutual_angle = random::Uniform(-20, 20) * 1_deg;

        auto d_orb = Elliptic(HD3167.mass, HD3167d.mass, 0.07757, 0.0, 89.3_deg + mutual_angle, 0.0, 0.0, isotherm);

        auto c_orb = Elliptic(HD3167.mass, HD3167c.mass, 0.1841, 0.05, 89.3_deg, 0.0, 0.0, isotherm);

        auto x_orb = Elliptic(HD3167.mass, HD3167x.mass, a_x, 0.0, 0.0, 0.0, 0.0, isotherm);

        move_particles(d_orb, HD3167d);
        move_particles(c_orb, HD3167c);
        move_particles(x_orb, HD3167x);

        move_to_COM_frame(HD3167, HD3167d, HD3167c, HD3167x);

        Solver sim{0, HD3167, HD3167d, HD3167c, HD3167x};

        Solver::RunArgs args;

        args.atol = 1e-13;

        double end_time_limit = 1e5_year;

        args.add_stop_condition(end_time_limit);

        int ejection_id = 0;

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
        args.add_stop_condition(any_ejection);
        // args.add_operation([&](auto& ptc, auto h) { std::cout << ptc.time() << '\n'; });
        args.add_stop_point_operation([&](auto& ptc, auto h) {
            print(file, ptc.time() / 1_year, ',', mass_x / 1_Mj, ',', a_x, ',', mutual_angle / 1_deg, ',', ejection_id,
                  '\n');
            file << std::endl;
        });

        sim.run(args);
    }
}

int main(int argc, char** argv) {
    size_t n = 100000;  // total scattering number
    size_t job_num = multi_thread::machine_thread_num;
    tf::Executor executor;

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num);
    }
    executor.wait_for_all();

    return 0;
}
