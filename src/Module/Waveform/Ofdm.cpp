#include "Module/Waveform/Ofdm.hpp"
#include <fftw3.h>
#include <algorithm>
#include <string.h>

using namespace aff3ct;
using namespace aff3ct::module;

template <typename B>
void Ofdm<B>::_fft(const B *X_K, B *Y_K)
{
    memcpy(fft_in, X_K, M * sizeof(std::complex<B>));
    fftwf_execute(fft_plan);
    memcpy(Y_K, fft_out, M * sizeof(std::complex<B>));
}

template <typename B>
void Ofdm<B>::_ifft(const B *X_K, B *Y_K)
{
    memcpy(ifft_in, X_K, M * sizeof(std::complex<B>));
    fftwf_execute(ifft_plan);
    memcpy(Y_K, ifft_out, M * sizeof(std::complex<B>));
}

#include "Tools/types.h"
template class aff3ct::module::Ofdm<R_32>;
