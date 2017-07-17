#include "Tools/Exception/exception.hpp"

#include "Module/Encoder/LDPC/Encoder_LDPC.hpp"
#include "Module/Encoder/LDPC/From_H/Encoder_LDPC_from_H.hpp"
#include "Module/Encoder/LDPC/DVBS2/Encoder_LDPC_DVBS2.hpp"

#include "Factory_encoder_LDPC.hpp"

using namespace aff3ct::module;
using namespace aff3ct::tools;

template <typename B>
Encoder_LDPC<B>* Factory_encoder_LDPC
::build(const parameters    &params,
        const Sparse_matrix &G,
        const Sparse_matrix &H)
{
	     if (params.type == "LDPC"      ) return new Encoder_LDPC       <B>(params.K, params.N_cw, G, params.n_frames);
	else if (params.type == "LDPC_H"    ) return new Encoder_LDPC_from_H<B>(params.K, params.N_cw, H, params.n_frames);
	else if (params.type == "LDPC_DVBS2") return new Encoder_LDPC_DVBS2 <B>(params.K, params.N_cw,    params.n_frames);

	throw cannot_allocate(__FILE__, __LINE__, __func__);
}

void Factory_encoder_LDPC
::build_args(Arguments_reader::arg_map &req_args, Arguments_reader::arg_map &opt_args)
{
	Factory_encoder::build_args(req_args, opt_args);

	opt_args[{"enc-type"}][2] += ", LDPC, LDPC_H, LDPC_DVBS2";

	opt_args[{"enc-h-path"}] =
		{"string",
		 "path to the H matrix (AList formated file, required by the \"LDPC_H\" encoder)."};

	opt_args[{"enc-g-path"}] =
		{"string",
		 "path to the G matrix (AList formated file, required by the \"LDPC\" encoder)."};
}

void Factory_encoder_LDPC
::store_args(const Arguments_reader& ar, parameters &params)
{
	params.type = "AZCW";

	Factory_encoder::store_args(ar, params);

	if(ar.exist_arg({"enc-h-path"})) params.H_alist_path = ar.get_arg({"enc-h-path"});
	if(ar.exist_arg({"enc-g-path"})) params.G_alist_path = ar.get_arg({"enc-g-path"});
}

void Factory_encoder_LDPC
::group_args(Arguments_reader::arg_grp& ar)
{
	Factory_encoder::group_args(ar);
}

void Factory_encoder_LDPC
::header(params_list& head_enc, const parameters& params)
{
	Factory_encoder::header(head_enc, params);

	if (params.type == "LDPC")
		head_enc.push_back(std::make_pair("G matrix path", params.G_alist_path));
	if (params.type == "LDPC_H")
		head_enc.push_back(std::make_pair("H matrix path", params.H_alist_path));
}

// ==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef MULTI_PREC
template aff3ct::module::Encoder_LDPC<B_8 >* aff3ct::tools::Factory_encoder_LDPC::build<B_8 >(const aff3ct::tools::Factory_encoder_LDPC::parameters&, const aff3ct::tools::Sparse_matrix&, const aff3ct::tools::Sparse_matrix&);
template aff3ct::module::Encoder_LDPC<B_16>* aff3ct::tools::Factory_encoder_LDPC::build<B_16>(const aff3ct::tools::Factory_encoder_LDPC::parameters&, const aff3ct::tools::Sparse_matrix&, const aff3ct::tools::Sparse_matrix&);
template aff3ct::module::Encoder_LDPC<B_32>* aff3ct::tools::Factory_encoder_LDPC::build<B_32>(const aff3ct::tools::Factory_encoder_LDPC::parameters&, const aff3ct::tools::Sparse_matrix&, const aff3ct::tools::Sparse_matrix&);
template aff3ct::module::Encoder_LDPC<B_64>* aff3ct::tools::Factory_encoder_LDPC::build<B_64>(const aff3ct::tools::Factory_encoder_LDPC::parameters&, const aff3ct::tools::Sparse_matrix&, const aff3ct::tools::Sparse_matrix&);
#else
template aff3ct::module::Encoder_LDPC<B>* aff3ct::tools::Factory_encoder_LDPC::build<B>(const aff3ct::tools::Factory_encoder_LDPC::parameters&, const aff3ct::tools::Sparse_matrix&, const aff3ct::tools::Sparse_matrix&);
#endif
// ==================================================================================== explicit template instantiation