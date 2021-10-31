#include "Module/Waveform/Ofdm.hpp"
#include <fftw3.h>
#include <algorithm>
#include <string.h>

using namespace aff3ct;
using namespace aff3ct::module;

template <typename B>
void Ofdm<B>::modulate(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K, int frame_id )
{
    for (size_t i = 0; i < N; i++)
    {
        memcpy(modu_in,X_K.data()+i*M,M*sizeof(std::complex<B>));
        fftw_execute(modu);
        for(size_t j=0;j<M;j++)
        {   
            Y_K.push_back(std::complex<B>(modu_out[j][0],modu_out[j][1]));
        }
    }
}

template <typename B>
void Ofdm<B>::demodulate(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K, int frame_id )
{
    for (size_t i = 0; i < N; i++)
    {
        memcpy(demodu_in,X_K.data()+i*M,M*sizeof(std::complex<B>));
        fftw_execute(demodu);
        for(size_t j=0;j<M;j++)
        {   
            Y_K.push_back(std::complex<B>(demodu_out[j][0],demodu_out[j][1]));
        }
    }
}

#include "Tools/types.h"
template class aff3ct::module::Ofdm<Q_64>;
