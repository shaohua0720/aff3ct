#include <iostream>

#include "../../Simulation/BFER/LDPC/Simulation_LDPC.hpp"
#include "../../Tools/bash_tools.h"

#include "Launcher_BFER_LDPC.hpp"

template <typename B, typename R, typename Q>
Launcher_BFER_LDPC<B,R,Q>
::Launcher_BFER_LDPC(const int argc, const char **argv, std::ostream &stream)
: Launcher_BFER<B,R,Q>(argc, argv, stream)
{
	// override parameters
	this->chan_params.quant_n_bits    = 6;
	this->chan_params.quant_point_pos = 2;

	// default parameters
	this->deco_params.max_iter        = 6;
	this->code_params.type            = "LDPC";
	this->deco_params.algo            = "BP_MIN_SUM";
	this->deco_params.implem          = "NAIVE";
}

template <typename B, typename R, typename Q>
void Launcher_BFER_LDPC<B,R,Q>
::build_args()
{
	Launcher_BFER<B,R,Q>::build_args();

	this->opt_args["max-iter"] = "n_iterations";
	this->doc_args["max-iter"] = "maximal number of iterations in the turbo decoder.";
}

template <typename B, typename R, typename Q>
void Launcher_BFER_LDPC<B,R,Q>
::store_args()
{
	Launcher_BFER<B,R,Q>::store_args();

	if(this->ar.exist_arg("max-iter")) this->deco_params.max_iter = std::stoi(this->ar.get_arg("max-iter"));
}

template <typename B, typename R, typename Q>
void Launcher_BFER_LDPC<B,R,Q>
::print_header()
{
	Launcher_BFER<B,R,Q>::print_header();

	this->stream << "# " << bold("* Decoding iterations per frame ") << " = " << this->deco_params.max_iter << std::endl;
}

template <typename B, typename R, typename Q>
void Launcher_BFER_LDPC<B,R,Q>
::build_simu()
{
	this->simu = new Simulation_LDPC<B,R,Q>(this->simu_params, 
	                                        this->code_params, 
	                                        this->enco_params, 
	                                        this->mod_params,
	                                        this->chan_params,
	                                        this->deco_params);
}

// ==================================================================================== explicit template instantiation 
#include "../../Tools/types.h"
#ifdef MULTI_PREC
template class Launcher_BFER_LDPC<B_8,R_8,Q_8>;
template class Launcher_BFER_LDPC<B_16,R_16,Q_16>;
template class Launcher_BFER_LDPC<B_32,R_32,Q_32>;
template class Launcher_BFER_LDPC<B_64,R_64,Q_64>;
#else
template class Launcher_BFER_LDPC<B,R,Q>;
#endif
// ==================================================================================== explicit template instantiation
