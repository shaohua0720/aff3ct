#ifndef DECODER_LDPC_BP_MIN_SUM_NAIVE_HPP_
#define DECODER_LDPC_BP_MIN_SUM_NAIVE_HPP_

#include "Decoder_LDPC_BP_naive.hpp"

template <typename B, typename R>
class Decoder_LDPC_BP_min_sum_naive : public Decoder_LDPC_BP_naive<B,R>
{
public:
	Decoder_LDPC_BP_min_sum_naive(const int &N, const int &K, const int& n_ite, 
	                              const std::vector<unsigned char> &n_variables_per_parity,
	                              const std::vector<unsigned char> &n_parities_per_variable,
	                              const std::vector<unsigned int > &transpose);
	virtual ~Decoder_LDPC_BP_min_sum_naive();

protected:
	// BP functions for decoding
	virtual bool BP_decode(const mipp::vector<R> &Y_N);
};

#endif /* DECODER_LDPC_BP_MIN_SUM_NAIVE_HPP_ */