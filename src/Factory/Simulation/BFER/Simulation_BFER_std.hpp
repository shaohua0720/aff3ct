#ifndef FACTORY_SIMULATION_BFER_STD_HPP_
#define FACTORY_SIMULATION_BFER_STD_HPP_

#include <string>

#include "Tools/Arguments_reader.hpp"
#include "Tools/Codec/Codec.hpp"

#include "Simulation_BFER.hpp"

namespace aff3ct
{
namespace simulation
{
template <typename B, typename R, typename Q>
class Simulation_BFER_std;
}
}

namespace aff3ct
{
namespace factory
{
struct Simulation_BFER_std : Simulation_BFER
{
	static const std::string name;
	static const std::string prefix;

	struct parameters : Simulation_BFER::parameters
	{
		parameters() : Simulation_BFER::parameters() {}

		virtual ~parameters() {}

		bool debug_fe = false;
	};

	template <typename B = int, typename R = float, typename Q = R>
	static simulation::Simulation_BFER_std<B,R,Q>* build(const parameters &params, tools::Codec<B,Q> &codec);

	static void build_args(arg_map &req_args, arg_map &opt_args, const std::string p = prefix);
	static void store_args(const tools::Arguments_reader& ar, parameters &params, const std::string p = prefix);
	static void header(params_list& head_sim, const parameters& params);
};
}
}

#endif /* FACTORY_SIMULATION_BFER_STD_HPP_ */