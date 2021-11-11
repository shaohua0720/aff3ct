#include <algorithm>
#include <sstream>
#include <string>
#include <cmath>
#include <complex>
#include <assert.h>

#include "Tools/Noise/Noise.hpp"
#include "Tools/Algo/Draw_generator/Gaussian_noise_generator/Standard/Gaussian_noise_generator_std.hpp"
#include "Tools/Algo/Draw_generator/Gaussian_noise_generator/Fast/Gaussian_noise_generator_fast.hpp"
#include "Tools/Exception/exception.hpp"
#include "Module/Channel/Raytrace/Raytrace_Generic.hpp"
#include "Module/Channel/Channel.hpp"

using namespace aff3ct;
using namespace aff3ct::module;

template <typename R>
float Raytrace_Generic<R>::initialTime = 0;

template <typename R>
tools::Gaussian_gen<R> *create_gaussian_generator(const tools::Gaussian_noise_generator_implem implem, const int seed)
{
    switch (implem)
    {
    case tools::Gaussian_noise_generator_implem::STD:
        return new tools::Gaussian_noise_generator_std<R>(seed);
        break;
    case tools::Gaussian_noise_generator_implem::FAST:
        return new tools::Gaussian_noise_generator_fast<R>(seed);
        break;
#ifdef AFF3CT_CHANNEL_GSL
    case tools::Gaussian_noise_generator_implem::GSL:
        return new tools::Gaussian_noise_generator_GSL<R>(seed);
        break;
#endif
#ifdef AFF3CT_CHANNEL_MKL
    case tools::Gaussian_noise_generator_implem::MKL:
        return new tools::Gaussian_noise_generator_MKL<R>(seed);
        break;
#endif
    default:
        std::stringstream message;
        message << "Unsupported 'implem' ('implem' = " << (int)implem << ").";
        throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
    };
}

template <typename R>
Raytrace_Generic<R>::Raytrace_Generic(int N, channel_model model, int sample_rate, float RMSDelaySpread,float dopplerHz,
                                      const tools::Gaussian_noise_generator_implem implem,
                                      int seed, bool gainNormalized) : Channel<R>(N), model(model), sample_rate(sample_rate),dopplerHz(dopplerHz),
                                                                       RMSDelaySpread(RMSDelaySpread), gainNormalized(gainNormalized),
                                                                       gaussian_generator(create_gaussian_generator<R>(implem, seed))
{
    this->initDelayGain();
}

template <typename R>
void Raytrace_Generic<R>::initDelayGain()
{
    switch (this->model)
    {
    case channel_model::EPA:
        this->delayProfile = {0, 30e-9, 70e-9, 90e-9, 110e-9, 190e-9, 410e-9};
        this->pathGaindB = {0.0, -1.0, -2.0, -3.0, -8.0, -17.2, -20.8};
        assert(this->delayProfile.size() == this->pathGaindB.size());
        this->numPath = this->delayProfile.size();
        break;

    case channel_model::EVA:
        this->delayProfile = {0, 30e-9, 150e-9, 310e-9, 370e-9, 710e-9, 1090e-9, 1730e-9, 2510e-9};
        this->pathGaindB = {0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9};
        assert(this->delayProfile.size() == this->pathGaindB.size());
        this->numPath = this->delayProfile.size();
        break;

    case channel_model::ETU:
        this->delayProfile = {0, 50e-9, 120e-9, 200e-9, 230e-9, 500e-9, 1600e-9, 2300e-9, 5000e-9};
        this->pathGaindB = {-1.0, -1.0, -1.0, 0.0, 0.0, 0.0, -3.0, -5.0, -7.0};
        assert(this->delayProfile.size() == this->pathGaindB.size());
        this->numPath = this->delayProfile.size();
        break;

    case channel_model::TDLA:
        this->delayProfile = {0, 0.3819, 0.4025, 0.5868, 0.461, 0.5375, 0.6708, 0.575, 0.7618, 1.5375, 1.8978, 2.2242,
                              2.1718, 2.4942, 2.5119, 3.0582, 4.081, 4.4579, 4.5695, 4.7966, 5.0066, 5.3043, 9.6586};
        this->pathGaindB = {-13.4, 0, -2.2, -4, -6, -8.2, -9.9, -10.5, -7.5, -15.9, -6.6, -16.7, -12.4, -15.2, -10.8,
                            -11.3, -12.7, -16.2, -18.3, -18.9, -16.6, -19.9, -29.7};
        assert(this->delayProfile.size() == this->pathGaindB.size());
        this->numPath = this->delayProfile.size();
        break;

    case channel_model::TDLB:
        this->delayProfile = {0, 0.1072, 0.2155, 0.2095, 0.287, 0.2986, 0.3752, 0.5055, 0.3681, 0.3697, 0.57, 0.5283,
                              1.1021, 1.2756, 1.5474, 1.7842, 2.0169, 2.8294, 3.0219, 3.6187, 4.1067, 4.279, 4.7834};
        this->pathGaindB = {0, -2.2, -4, -3.2, -9.8, -1.2, -3.4, -5.2, -7.6, -3, -8.9, -9, -4.8, -5.7, -7.5, -1.9, -7.6,
                            -12.2, -9.8, -11.4, -14.9, -9.2, -11.3};
        assert(this->delayProfile.size() == this->pathGaindB.size());
        this->numPath = this->delayProfile.size();
        break;

    case channel_model::TDLC:
        this->delayProfile = {0, 0.2099, 0.2219, 0.2329, 0.2176, 0.6366, 0.6448, 0.656, 0.6584, 0.7935,
                              0.8213, 0.9336, 1.2285, 1.3083, 2.1704, 2.7105, 4.2589, 4.6003, 5.4902, 5.6077, 6.3065, 6.6374, 7.0427, 8.6523};
        this->pathGaindB = {-4.4, -1.2, -3.5, -5.2, -2.5, 0, -2.2, -3.9, -7.4, -7.1, -10.7, -11.1, -5.1, -6.8, -8.7,
                            -13.2, -13.9, -13.9, -15.8, -17.1, -16, -15.7, -21.6, -22.8};
        assert(this->delayProfile.size() == this->pathGaindB.size());
        this->numPath = this->delayProfile.size();
        break;

    case channel_model::User:
        this->delayProfile = {};
        this->pathGaindB = {};
        this->numPath = this->delayProfile.size();
        break;
    default:
        break;
    }
    if ((this->model == channel_model::TDLA) || (this->model == channel_model::TDLB) || (this->model == channel_model::TDLC))
    {
        for (auto &i : this->delayProfile)
        {
            i = i * this->RMSDelaySpread * 1e-9;
        }
    }

    if (gainNormalized)
    {
        auto power = 0.0;
        for (auto i : this->pathGaindB)
        {
            power += std::pow(10, i / 10.0);
        }
        for (auto &i : this->pathGaindB)
        {
            auto t = std::pow(10, i / 10.0) / power;
            this->pathGain.push_back(t);
            i = 10 * log10(t);
        }
    }

    for(int i =0;i<this->numPath;i++)
    {
        std::random_device rd;
        std::default_random_engine rd_engine(rd());
        std::uniform_real_distribution<float> angle_t(0,1);
        std::normal_distribution<float> randn(0, 1);

        auto angle = angle_t(rd_engine)*2.0f*M_PI;
        //std::cout << "################" << angle << std::endl;
        this->AoA.push_back(angle);
        this->pathDoppler.push_back(this->dopplerHz*cos(angle));
        this->pathOffset.push_back(2*M_PI*angle_t(rd_engine));

        std::complex<double> gain = sqrt(this->pathGain[i])/sqrt(2.0)*std::complex<double>(randn(rd_engine),randn(rd_engine));
        //std::cout<<"############## complex gain:"<<gain<<" "<<std::endl;
        this->pathRayleighGain.push_back(gain);
    }

    for(auto& i : this->delayProfile)
    {
        this->delayInSample.push_back(ceil(i*this->sample_rate)); //floor or ceil or round
        i = this->delayInSample.back()*1.0d/this->sample_rate;
    }
}

template <typename R>
void Raytrace_Generic<R>::_add_noise(const float *CP, const R *X_N, R *Y_N, const size_t frame_id)
{
    //std::cout << "-------- num of "<<this->N<<" real number in total"<<std::endl;
    std::vector<R> awgn_noise(this->N);
    gaussian_generator->generate(awgn_noise.data(), this->N, (R)*CP);

    for(size_t i =0;i<this->N;)
    {
        std::complex<R> out_point(0,0);
        std::complex<R> in_point = std::complex<R>(X_N[i],X_N[i+1]);
        for(auto j =0;j< this->numPath;j++)
        {
            double time = (i/2.0-this->delayInSample[j])/this->sample_rate+this->initialTime;
            auto phase = exp(std::complex<R>(0,2*M_PI*(this->pathDoppler[j]*time+this->pathOffset[j])));
            out_point += (std::complex<R>)this->pathRayleighGain[j]*phase*in_point;
        }
        // add awgn noise
        out_point += std::complex<R>(1.0f/sqrt(2)*awgn_noise[i],1.0f/sqrt(2)*awgn_noise[i+1]);

        Y_N[i]=out_point.real();
        Y_N[i+1]=out_point.imag();
        i+=2;
    }
    this->initialTime+=this->N/2;
}

#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Raytrace_Generic<R_32>;
#else
template class aff3ct::module::Raytrace_Generic<R_32>;
#endif