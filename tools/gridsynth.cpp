#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

#include "nwqec/gridsynth/gridsynth.hpp"
#include "nwqec/core/constants.hpp"

void print_usage(const char *program_name)
{
    std::cout << "Usage: " << program_name << " <angle> [epsilon]" << std::endl;
    std::cout << "  angle   - Target rotation angle in radians (or 'pi/n')" << std::endl;
    std::cout << "  epsilon - Optional absolute tolerance (e.g., 1e-6).\n"
                 "            If omitted, defaults to |theta| * 1e-2." << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  " << program_name << " 0.785398 1e-10   # π/4 with ε=1e-10" << std::endl;
    std::cout << "  " << program_name << " pi/4             # ε defaults to |θ|*1e-2" << std::endl;
    std::cout << "  " << program_name << " pi/8 1e-12       # π/8 with ε=1e-12" << std::endl;
}

int main(int argc, char *argv[])
{
    std::string theta;
    std::string epsilon;

    if (argc < 2)
    {
        print_usage(argv[0]);
        return 1;
    }

    theta = argv[1];

    if (argc >= 3)
    {
        // Use provided epsilon verbatim (supports values like 1e-10 or 0.001)
        epsilon = argv[2];
    }
    else
    {
        // Compute epsilon = |theta| * 1e-2
        double theta_abs = 0.0;
        try
        {
            // Try numeric parse first
            theta_abs = std::fabs(std::stod(theta));
        }
        catch (...)
        {
            // Fallback: parse simple 'pi/<n>' forms
            const std::string prefix = "pi/";
            auto pos = theta.find(prefix);
            if (pos != std::string::npos)
            {
                std::string denom_str = theta.substr(pos + prefix.size());
                try
                {
                    double denom = std::stod(denom_str);
                    if (denom != 0.0)
                    {
                        theta_abs = M_PI / std::fabs(denom);
                    }
                }
                catch (...)
                {
                    theta_abs = 0.0;
                }
            }
        }
        double eps_val = (theta_abs > 0.0) ? theta_abs * NWQEC::DEFAULT_EPSILON_MULTIPLIER : 1e-10;
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(10) << eps_val;
        epsilon = oss.str();
    }

    // Print input parameters
    std::cout << "Gridsynth Parameters:" << std::endl;
    std::cout << "  θ (theta) = " << theta
              << " rad"
              << std::endl;

    std::cout << "  ε (epsilon) = " << epsilon << std::endl;
    std::cout << std::endl;

    // Call gridsynth
    std::cout << "Computing optimal gate sequence..." << std::endl;

    // Use high precision values for computation
    auto gates = gridsynth::gridsynth_gates(theta, epsilon,
                                            NWQEC::DEFAULT_DIOPHANTINE_TIMEOUT_MS,
                                            NWQEC::DEFAULT_FACTORING_TIMEOUT_MS,
                                            false, true);

    // Print results
    std::cout << std::endl;
    std::cout << "Results:" << std::endl;
    std::cout << "  Gate sequence: " << gates << std::endl
              << std::endl;
    std::cout << "  Gate count: " << gates.size() << std::endl;

    // Calculate and display actual error
    std::string actual_error = gridsynth::error(theta, gates);
    std::cout << "  Actual error: " << actual_error << std::endl;
    std::cout << "  Target error:  " << epsilon << std::endl;

    // Count gate types
    int t_count = 0, h_count = 0, s_count = 0, w_count = 0;
    for (const auto &gate : gates)
    {
        if (gate == 'T')
            t_count++;
        else if (gate == 'H')
            h_count++;
        else if (gate == 'S')
            s_count++;
        else if (gate == 'W')
            w_count++;
    }

    std::cout << std::endl;
    std::cout << "Gate breakdown:" << std::endl;
    if (t_count > 0)
        std::cout << "  T gates: " << t_count << std::endl;
    if (h_count > 0)
        std::cout << "  H gates: " << h_count << std::endl;
    if (s_count > 0)
        std::cout << "  S gates: " << s_count << std::endl;
    if (w_count > 0)
        std::cout << "  W gates: " << w_count << std::endl;

    return 0;
}
