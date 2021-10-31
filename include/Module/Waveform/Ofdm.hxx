#include "Tools/Exception/exception.hpp"
#include "Module/Waveform/Ofdm.hpp"
#include <numeric>
#include <math.h>

namespace aff3ct
{
    namespace module
    {
        template <typename B>
        Task &Ofdm<B>::operator[](const waveform::tsk t)
        {
            return Module::operator[]((size_t)t);
        }

        template <typename B>
        Socket &Ofdm<B>::operator[](const waveform::sck::modulate s)
        {
            return Module::operator[]((size_t)waveform::tsk::modulate)[(size_t)s];
        }

        template <typename B>
        Socket &Ofdm<B>::operator[](const waveform::sck::demodulate s)
        {
            return Module::operator[]((size_t)waveform::tsk::modulate)[(size_t)s];
        }

        template <typename B>
        Ofdm<B>::Ofdm(const int M, const int N, const bool padding) : Module(), M(M), N(N), padding(padding)
        {

            const std::string name = "OFDM";
            this->set_name(name);
            this->set_short_name(name);

            if (padding)
            {
                fft_size = pow(2, ceil(log2(M)));
                start_pos = floor((fft_size - M) / 2); // place at the center band.
                end_pos = start_pos + M - 1;
            }
            else
            {
                fft_size = M;
            }

            std::cout << "fft size:" << M << std::endl;

            fft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);
            fft_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);

            ifft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);
            ifft_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);

            fft_plan = fftw_plan_dft_1d(fft_size, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE);
            ifft_plan = fftw_plan_dft_1d(fft_size, ifft_in, ifft_out, FFTW_BACKWARD, FFTW_MEASURE);

            // create related modulate and demodulates
            auto &p_mod = this->create_task("modulate");
            auto pin_X_K = this->template create_socket_in<B>(p_mod, "X_K", 2 * M * N); //std::complex<B> = 2*B
            auto pin_Y_K = this->template create_socket_out<B>(p_mod, "Y_K", 2 * M * N);

            this->create_codelet(p_mod, [pin_X_K, pin_Y_K](Module &m, Task &t, const size_t frame_id) -> int
                                 {
                                     auto &ofdm = static_cast<Ofdm<B> &>(m);
                                     B *in = static_cast<B *>(t[pin_X_K].get_dataptr());
                                     B *out = static_cast<B *>(t[pin_Y_K].get_dataptr());
                                     return status_t::SUCCESS;
                                 });
        }

        template <typename B>
        void Ofdm<B>::setCPlength(std::vector<int> cp)
        {
            assert(cp.size() == N);
            this->cp.insert(this->cp.begin(), cp.begin(), cp.end());
        }

        template <typename B>
        Ofdm<B>::~Ofdm()
        {
            fftw_destroy_plan(fft_plan); // destroy the fftw
            fftw_destroy_plan(ifft_plan);
            fftw_free(fft_in);
            fftw_free(fft_out);
            fftw_free(ifft_in);
            fftw_free(ifft_out);
        }

        template <typename B>
        Ofdm<B> *Ofdm<B>::clone() const
        {
            throw tools::unimplemented_error(__FILE__, __LINE__, __func__);
        }

        template <typename B>
        void Ofdm<B>::modulate(const B *X_K, B *Y_K, int frame_id)
        {
            B data[fft_size * sizeof(std::complex<B>)]; // for one symbol
            size_t pos = 0;
            for (size_t i = 0; i < N; i++)
            {
                this->_ifft(X_K + i * M * sizeof(std::complex<B>), Y_K + pos + cp[i]);
                memcpy(Y_K + pos, Y_K + pos + fft_size, cp[i] * sizeof(std::complex<B>)); //copy the cp
                pos = pos + cp[i] + fft_size;                                             //update the cursor
            }
        }

        template <typename B>
        void Ofdm<B>::demodulate(const B *X_K, B *Y_K, int frame_id)
        {
            size_t pos = 0;
            for (size_t i = 0; i < N; i++)
            {
                this->_fft(X_K + pos + cp[i], Y_K + i * M * sizeof(std::complex<B>)); //discard cp
                pos = pos + cp[i] + fft_size;                                         //update the cursor
            }
        }
    }
}