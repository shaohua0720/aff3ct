#include "Module/Waveform/Ofdm.hpp"
#include <fftw3.h>
#include <algorithm>

using namespace aff3ct;
using namespace aff3ct::module;

template <typename B>
void Ofdm<B>::modulate(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K, int frame_id )
{
    for (size_t i = 0; i < N; i++)
    {
        std::copy_n(X_K.begin()+i*M,M,modu_in);
        fftw_execute(modu);
        std::copy_n(modu_out,M,Y_K.begin()+i*M);
    }
}

template <typename B>
void Ofdm<B>::demodulate(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K, int frame_id )
{
    for (size_t i = 0; i < N; i++)
    {
        std::copy_n(X_K.begin()+i*M,M,demodu_in);
        fftw_execute(modu);
        std::copy_n(demodu_out,M,Y_K.begin()+i*M);
    }
}
