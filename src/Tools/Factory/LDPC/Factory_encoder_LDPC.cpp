#include "../../../Encoder/LDPC/Encoder_LDPC_fake.hpp"

#include "Factory_encoder_LDPC.hpp"

template <typename B>
Encoder_LDPC_sys<B>* Factory_encoder_LDPC<B>
::build(const t_simulation_param &simu_params,
        const t_code_param       &code_params,
        const t_encoder_param    &enco_params)
{
	Encoder_LDPC_sys<B> *encoder = nullptr;

	// build the encoder
	if (enco_params.systematic)
		encoder = new Encoder_LDPC_fake<B>(code_params.K, code_params.N);

	return encoder;
}

// ==================================================================================== explicit template instantiation 
#include "../../types.h"
#ifdef MULTI_PREC
template struct Factory_encoder_LDPC<B_8>;
template struct Factory_encoder_LDPC<B_16>;
template struct Factory_encoder_LDPC<B_32>;
template struct Factory_encoder_LDPC<B_64>;
#else
template struct Factory_encoder_LDPC<B>;
#endif
// ==================================================================================== explicit template instantiation