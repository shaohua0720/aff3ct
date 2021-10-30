/**
 * 
 * 
 * 
 **/
#ifndef OFDM_HPP_
#define OFDM_HPP_

#include <iostream>
#include <complex>
#include <fftw3.h>
#include "Tools/Interface/Interface_set_seed.hpp"
#include "Tools/Interface/Interface_is_done.hpp"
#include "Tools/Interface/Interface_reset.hpp"
#include "Module/Task.hpp"
#include "Module/Socket.hpp"
#include "Module/Module.hpp"

namespace aff3ct
{
    namespace module
    {
        namespace waveform
        {
            enum class tsk : size_t
            {
                modulate, demodulate, SIZE
            };
            namespace sck
            {
                enum class modulate : size_t
                {
                    X_K,
                    status
                };
                enum class demodulate : size_t
                {
                    Y_K,
                    status
                };
            }
        }

        template <typename B = double>
        class Ofdm : public Module
        {
        public:
            inline Task &operator[](const waveform::tsk t);
            inline Socket &operator[](const waveform::sck::modulate s);
            inline Socket &operator[](const waveform::sck::demodulate s);

        protected:
            const int M, N; /* Number of frequency and time domain bins. */

        public:
            Ofdm(const int M, const int N = 14, const bool padding=true);
            virtual ~Ofdm();
            virtual Ofdm<B> *clone() const;

            void setCPlength(const std::vector<int> cp);

            void modulate(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K, int frame_id = -1);
            void demodulate(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K, int frame_id = -1);
        private:
            std::vector<int> cp; // cyclic prefix points
            int fft_size;
            bool padding; //if padding is enabled for fft and scs numbers.
            int start_pos=0;
            int end_pos=M-1;

            fftw_plan modu;
            fftw_complex *modu_in, *modu_out;
            fftw_plan demodu;
            fftw_complex *demodu_in, *demodu_out;
        };
    }

} // namespace aff3ct

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Module/Waveform/Ofdm.hxx"
#endif

#endif