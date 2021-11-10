#include "Module/Waveform/OTFS/Otfs.hpp"
#include <complex>
#include <string>
#include <iostream>

namespace aff3ct
{
    namespace module
    {
        template <typename B>
        Otfs<B>::Otfs(int M, const std::vector<int> cp, int N):Ofdm<B>(M,cp,N)
        {
            const std::string name = "OTFS";
            this->set_name(name);

            for(auto i : cp)
                cp_cplx += i;

            isfft_data_inplace = (fftwf_complex *)fftwf_malloc(M*N*sizeof(fftwf_complex));
            /*N-points IFFT for doppler*/
            isfft_doppler = fftwf_plan_many_dft(1,&N,M,isfft_data_inplace,NULL,M,1,isfft_data_inplace,NULL,M,1,FFTW_BACKWARD,FFTW_MEASURE);
            /*M-points FFT for delay*/
            isfft_delay = fftwf_plan_many_dft(1,&M,N,isfft_data_inplace,NULL,1,M,isfft_data_inplace,NULL,1,M,FFTW_FORWARD,FFTW_MEASURE);
            
            sfft_data_inplace = (fftwf_complex *)fftwf_malloc(M*N*sizeof(fftwf_complex));
            /*N-points FFT for doppler*/
            sfft_doppler = fftwf_plan_many_dft(1,&N,M,sfft_data_inplace,NULL,M,1,sfft_data_inplace,NULL,M,1,FFTW_FORWARD,FFTW_MEASURE);
            /*M-points IFFT for delay*/
            sfft_delay = fftwf_plan_many_dft(1,&M,N,sfft_data_inplace,NULL,1,M,sfft_data_inplace,NULL,1,M,FFTW_BACKWARD,FFTW_MEASURE);

        }
        
        template <typename B>
        Otfs<B>::~Otfs()
        {
            fftwf_destroy_plan(isfft_delay);
            fftwf_destroy_plan(isfft_doppler);
            fftwf_destroy_plan(sfft_delay);
            fftwf_destroy_plan(sfft_doppler);
            fftwf_free(isfft_data_inplace);
            fftwf_free(sfft_data_inplace);
        }
        template <typename B>
        void Otfs<B>::modulate(const B *X_K, B *Y_K,int frame_id)
        {
            size_t count = this->M * this->N * sizeof(B) * 2;
            B* data = (B*)malloc(count);
            this->_isfft(X_K,data); // ISFFT transform
            this->_modulate(data,Y_K,frame_id); // Heisenberg transform
            this->normalize(Y_K, count/sizeof(B)+cp_cplx*2,1.0/(this->M*this->N)*sqrt(this->N));
            free(data);
        }
        template <typename B>
        void Otfs<B>::demodulate(const B *X_K, B *Y_K,int frame_id)
        {
            size_t count = this->M * this->N * sizeof(B) * 2;
            B* data = (B *)malloc(count);
            this->_demodulate(X_K,data,frame_id);  //Wiener transform 
            this->_sfft(data,Y_K);  // SFFT transform
            this->normalize(Y_K, count/sizeof(B),1.0/(this->M*sqrt(this->N)));
            free(data);
        }

        template <typename B>
        void Otfs<B>::_isfft(const B *X_K, B *Y_K) // IFFT for doppler(row), FFT for delay(col)
        {
            memcpy(isfft_data_inplace,X_K,this->M*this->N*sizeof(B)*2);
            fftwf_execute(isfft_delay);
            fftwf_execute(isfft_doppler);
            memcpy(Y_K,isfft_data_inplace,this->M*this->N*sizeof(B)*2);
        }
        template <typename B>
        void Otfs<B>::_sfft(const B *X_K, B *Y_K) // FFT for doppler(row), IFFT for delay(col)
        {  
            memcpy(sfft_data_inplace,X_K,this->M*this->N*sizeof(B)*2);
            fftwf_execute(sfft_doppler);
            fftwf_execute(sfft_delay);
            memcpy(Y_K,sfft_data_inplace,this->M*this->N*sizeof(B)*2);
        }

        template <typename B>
        void Otfs<B>::normalize(B* Y_K, size_t count,float weight)
        {
            for (auto i = 0; i < count; i++)
            {
                //std::cout<<i<<":"<<Y_K[i]<<std::endl;
                Y_K[i] = (B)(Y_K[i]*weight);
            }
        }
    }
}

#include "Tools/types.h"
template class aff3ct::module::Otfs<R_32>;