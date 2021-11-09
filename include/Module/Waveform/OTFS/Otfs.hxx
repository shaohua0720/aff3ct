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

            isfft_data_inplace = (fftw_complex *)fftw_malloc(M*N*sizeof(fftw_complex));
            /*N-points IFFT for doppler*/
            isfft_doppler = fftw_plan_many_dft(1,&N,M,isfft_data_inplace,NULL,M,1,isfft_data_inplace,NULL,M,1,FFTW_BACKWARD,FFTW_MEASURE);
            /*M-points FFT for delay*/
            isfft_delay = fftw_plan_many_dft(1,&M,N,isfft_data_inplace,NULL,1,M,isfft_data_inplace,NULL,1,M,FFTW_FORWARD,FFTW_MEASURE);
            
            sfft_data_inplace = (fftw_complex *)fftw_malloc(M*N*sizeof(fftw_complex));
            /*N-points FFT for doppler*/
            sfft_doppler = fftw_plan_many_dft(1,&N,M,sfft_data_inplace,NULL,M,1,sfft_data_inplace,NULL,M,1,FFTW_FORWARD,FFTW_MEASURE);
            /*M-points IFFT for delay*/
            sfft_delay = fftw_plan_many_dft(1,&M,N,sfft_data_inplace,NULL,1,M,sfft_data_inplace,NULL,1,M,FFTW_BACKWARD,FFTW_MEASURE);

        }
        
        template <typename B>
        Otfs<B>::~Otfs()
        {
            fftw_destroy_plan(isfft_delay);
            fftw_destroy_plan(isfft_doppler);
            fftw_destroy_plan(sfft_delay);
            fftw_destroy_plan(sfft_doppler);
            fftw_free(isfft_data_inplace);
            fftw_free(sfft_data_inplace);
        }
        template <typename B>
        void Otfs<B>::modulate(const B *X_K, B *Y_K,int frame_id)
        {
            size_t count = this->M*this->N*sizeof(fftw_complex);
            B* data = (B*)malloc(count * sizeof(B));
            //B data[count];
            this->_isfft(X_K,data); // ISFFT transform
            this->display_arrary(data,count/sizeof(B),"ISFFT:");
            this->_modulate(data,Y_K,frame_id); // Heisenberg transform
            free(data);
        }
        template <typename B>
        void Otfs<B>::demodulate(const B *X_K, B *Y_K,int frame_id)
        {
            B* data = (B*)malloc(this->M * this->N * sizeof(fftw_complex));
            //B data[this->M*this->N*sizeof(fftw_complex)];
            this->_demodulate(X_K,data,frame_id);  //wiener transform 
            this->_sfft(data,Y_K);  // SFFT transform
            free(data);
        }

        template <typename B>
        void Otfs<B>::_isfft(const B *X_K, B *Y_K) // IFFT for doppler(row), FFT for delay(col)
        {
            memcpy(isfft_data_inplace,X_K,this->M*this->N*sizeof(fftw_complex));
            fftw_execute(isfft_doppler);
            fftw_execute(isfft_delay);
            memcpy(Y_K,isfft_data_inplace,this->M*this->N*sizeof(fftw_complex));
        }
        template <typename B>
        void Otfs<B>::_sfft(const B *X_K, B *Y_K) // FFT for doppler(row), IFFT for delay(col)
        {  
            memcpy(sfft_data_inplace,X_K,this->M*this->N*sizeof(fftw_complex));
            fftw_execute(sfft_doppler);
            fftw_execute(sfft_delay);
            memcpy(Y_K,sfft_data_inplace,this->M*this->N*sizeof(fftw_complex));
        }
    }
}

#include "Tools/types.h"
template class aff3ct::module::Otfs<Q_64>;