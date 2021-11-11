#ifndef RAYTRACE_GENERIC_HPP_
#define RAYTRACE_GENERIC_HPP_

#include <memory>
#include <vector>
#include <complex>

#include "Tools/Algo/Draw_generator/Gaussian_noise_generator/Gaussian_noise_generator.hpp"
#include "Module/Channel/Channel.hpp"

namespace aff3ct
{
    namespace module
    {
        enum class channel_model : size_t
        {
            EPA,
            EVA,
            ETU,
            TDLA,
            TDLB,
            TDLC,
            User
        };

        template <typename R>
        class Raytrace_Generic : public Channel<R>
        {
        public:
            Raytrace_Generic(int N, channel_model model, int sample_rate, float RMSDelaySpread, float dopplerHz,
                             const tools::Gaussian_noise_generator_implem implem = tools::Gaussian_noise_generator_implem::STD,
                             int seed = 0, bool gainNormalized = true);
            
        protected:
            void _add_noise(const float *CP, const R *X_N, R *Y_N, const size_t frame_id);

        private:
            long int sample_rate;
            float RMSDelaySpread;
            const float dopplerHz;
            channel_model model;
            std::vector<double> delayProfile;
            std::vector<int>    delayInSample;
            std::vector<double> pathGaindB;
            std::vector<double> pathGain; // in non-octa
            std::vector<double> AoA;
            std::vector<double> pathDoppler;
            std::vector<std::complex<double>> pathRayleighGain; 
            std::vector<double> pathOffset; //phase offset for each path
            int numPath;
            bool gainNormalized;
            std::shared_ptr<tools::Gaussian_noise_generator<R>> gaussian_generator;
            static float initialTime; // for time evolution

        private:
            void initDelayGain();
        };
    }
}

#endif