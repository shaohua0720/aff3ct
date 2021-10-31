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

            const std::string name = "Waveform";
            this->set_name(name);
            this->set_short_name(name);

            if (padding)
            {
                fft_size = pow(2 , ceil(log2(M)));
                start_pos = floor((fft_size - M) / 2); // place at the center band.
                end_pos = start_pos + M - 1;
            }
            else
            {
                fft_size = M;
            }

            std::cout<<"fft size:"<<M<<std::endl;

            fft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);
            fft_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);

            ifft_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);
            ifft_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fft_size);

            fft_plan = fftw_plan_dft_1d(fft_size, fft_in, fft_out, FFTW_FORWARD, FFTW_MEASURE);
            ifft_plan = fftw_plan_dft_1d(fft_size, ifft_in, ifft_out, FFTW_BACKWARD, FFTW_MEASURE);

            // create related modulate and demodulates
            /* auto& p = this->create_task("modulate");
            auto ps_X_K = this->template create_socket_in<B>(p, "X_K", 1);
            this->create_codelet(p, [ps_X_K](Module &m, Task &t, const size_t frame_id) -> int
            {
                auto &ofdm = static_cast<Ofdm<B>&>(m);
            } */
        }

        template <typename B>
        void Ofdm<B>::setCPlength(std::vector<int> cp)
        {
            assert(cp.size()==N);
            this->cp.insert(this->cp.begin(),cp.begin(),cp.end());
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
        void Ofdm<B>::modulate(const std::vector<std::complex<B>> &constell, std::vector<std::complex<B>> &Y_K,int frame_id)
        {
            assert(constell.size()==M*N);
            static std::vector<std::complex<B>> mod_symWoCP;
            static std::vector<std::complex<B>> TDSym; // time domain points.
            
            std::cout<<std::endl;
            for(size_t i = 0 ;i < N; i++)
            {
                mod_symWoCP.insert(mod_symWoCP.begin()+start_pos,constell.begin()+i*M,constell.begin()+i*M+M);
                this->_ifft(mod_symWoCP,TDSym);
                
                Y_K.insert(Y_K.end(),TDSym.end()-this->cp[i],TDSym.end()); //insert CP
                Y_K.insert(Y_K.end(),TDSym.begin(),TDSym.end());
                mod_symWoCP.clear();
                TDSym.clear();
            }
        }

        template <typename B>
        void Ofdm<B>::demodulate(const std::vector<std::complex<B>> &X_K, std::vector<std::complex<B>> &Y_K,int frame_id)
        {
            static std::vector<std::complex<B>> demod_symWoCP;
            static std::vector<std::complex<B>> FDSym;
            auto idx = X_K.begin();
            for(size_t i = 0;i < N; i++)
            {
                demod_symWoCP.insert(demod_symWoCP.begin(),idx+this->cp[i],idx+this->cp[i]+fft_size); // discard cp
                idx = idx + this->cp[i]+fft_size;

                this->_fft(demod_symWoCP,FDSym);

                Y_K.insert(Y_K.end(),FDSym.begin()+start_pos,FDSym.begin()+start_pos+M); // extract signal
                demod_symWoCP.clear();
                FDSym.clear();
            }
        }
    }
}