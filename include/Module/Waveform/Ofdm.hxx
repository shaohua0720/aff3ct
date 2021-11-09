#include "Tools/Exception/exception.hpp"
#include "Module/Waveform/Ofdm.hpp"
#include <numeric>
#include <math.h>
#include <string.h>
#include <string>
#include <memory>

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
        Ofdm<B>::Ofdm(const int M, const std::vector<int> cps, const int N,const bool padding) : 
        Module(), 
        M(M), 
        N(N), 
        padding(padding),
        cp(cps)
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

            cp_sum = accumulate(cp.begin(), cp.end(), 0);

            fft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);
            fft_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);

            ifft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);
            ifft_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);

            fft_plan = fftw_plan_dft_1d(fft_size, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE);
            ifft_plan = fftw_plan_dft_1d(fft_size, ifft_in, ifft_out, FFTW_BACKWARD, FFTW_MEASURE);

            // create related modulate and demodulates
            auto &p_mod = this->create_task("modulate");
            auto pmod_X_K = this->template create_socket_in<B>(p_mod, "X_K", 2 * M * N); //std::complex<B> = 2*B
            auto pmod_Y_K = this->template create_socket_out<B>(p_mod, "Y_K", 2 * M * N+2*cp_sum);

            this->create_codelet(p_mod, [pmod_X_K, pmod_Y_K](Module &m, Task &t, const size_t frame_id) -> int
                                 {
                                     auto &ofdm = static_cast<Ofdm<B> &>(m);
                                     B *in = static_cast<B *>(t[pmod_X_K].get_dataptr());
                                     B *out = static_cast<B *>(t[pmod_Y_K].get_dataptr());
                                     ofdm.modulate(in,out);
                                     return status_t::SUCCESS;
                                 });

            auto &p_demod = this->create_task("demodulate");
            auto pdemod_X_K = this->template create_socket_in<B>(p_demod, "X_K", 2 * M * N+2*cp_sum);
            auto pdemod_Y_K = this->template create_socket_out<B>(p_demod, "Y_K", 2 * M * N);

            this->create_codelet(p_demod, [pdemod_X_K, pdemod_Y_K](Module &m, Task &t, const size_t frame_id) -> int
                                 {
                                     auto &ofdm = static_cast<Ofdm<B> &>(m);
                                     B *in = static_cast<B *>(t[pdemod_X_K].get_dataptr());
                                     B *out = static_cast<B *>(t[pdemod_Y_K].get_dataptr());
                                     ofdm.demodulate(in,out);
                                     return status_t::SUCCESS;
                                 });
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
        void Ofdm<B>::display_arrary(B *data, size_t size, const std::string info)
        {
            std::cout << std::endl
                      << info << ":" << std::endl;
            for (int i = 0; i < size; i++)
            {
                std::cout << data[i] << " ";
            }
            std::cout << std::endl;
        }

        template <typename B>
        void Ofdm<B>::modulate(const B *X_K, B *Y_K,int frame_id)
        {
            this->_modulate(X_K,Y_K,frame_id);
        }
        template <typename B>
        void Ofdm<B>::demodulate(const B *X_K,  B *Y_K,int frame_id)
        {
            this->_demodulate(X_K,Y_K,frame_id);
        }

        template <typename B>
        void Ofdm<B>::_modulate(const B *X_K, B *Y_K, int frame_id)
        {
            size_t B_w = sizeof(B);
            B* data = (B*)malloc(2 * fft_size * B_w);
            //B data[2 * fft_size];
            size_t pos = 0;
            for (size_t i = 0; i < N; i++)
            {
                memset(data, 0, fft_size * 2 * B_w);
                memcpy(data + 2 * start_pos * B_w, X_K + i * 2 * M, 2 * M * B_w);
                this->_ifft(data, Y_K + 2 * pos + 2 * cp[i]);
                memcpy(Y_K + pos * 2, Y_K + (pos + fft_size) * 2, cp[i] * 2 * B_w);
                pos = pos + cp[i] + fft_size; //update the cursor
            }
            free(data);
        }

        template <typename B>
        void Ofdm<B>::_demodulate(const B *X_K, B *Y_K, int frame_id)
        {
            size_t B_w = sizeof(B);
            B* data = (B*)malloc(2 * fft_size * B_w);
            //B data[2 * fft_size];
            size_t pos = 0;
            for (size_t i = 0; i < N; i++)
            {
                memset(data, 0, 2 * fft_size * B_w);
                this->_fft(X_K + (pos + cp[i]) * 2, data); //discard cp
                memcpy(Y_K + i * M * 2, data + start_pos * 2, 2*M * B_w);
                pos = pos + cp[i] + fft_size; //update the cursor
            }
            free(data);
        }
    }
}