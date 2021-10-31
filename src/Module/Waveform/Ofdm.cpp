#include "Module/Waveform/Ofdm.hpp"
#include <fftw3.h>
#include <algorithm>
#include <string.h>

using namespace aff3ct;
using namespace aff3ct::module;

template <typename B>
void Ofdm<B>::_fft(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K)
{
    memcpy(fft_in, X_K.data(), M * sizeof(std::complex<B>));
    fftw_execute(fft_plan);
    for (size_t j = 0; j < M; j++)
    {
        Y_K.push_back(std::complex<B>(fft_out[j][0], fft_out[j][1]));
    }
}

template <typename B>
void Ofdm<B>::_ifft(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K)
{
    memcpy(ifft_in, X_K.data(), M * sizeof(std::complex<B>));
    fftw_execute(ifft_plan);
    for (size_t j = 0; j < M; j++)
    {
        Y_K.push_back(std::complex<B>(ifft_out[j][0], ifft_out[j][1]));
    }
}

#include "Tools/types.h"
template class aff3ct::module::Ofdm<Q_64>;
