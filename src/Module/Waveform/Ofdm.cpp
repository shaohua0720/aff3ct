#include "Module/Waveform/Ofdm.hpp"
using namespace aff3ct;
using namespace aff3ct::module;

template <typename B>
template <class AB>
void Ofdm<B>::modulate(const std::vector<B, AB> &X_K, std::vector<B, AB> &Y_K, int frame_id)
{
    Y_K = X_K;
}

template <typename B>
template <class AB>
void Ofdm<B>::demodulate(const std::vector<B, AB> &X_K, std::vector<B, AB> &Y_K, int frame_id)
{
    Y_K = X_K;
}
