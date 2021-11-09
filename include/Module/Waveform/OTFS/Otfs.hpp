#ifndef OTFS_HPP_
#define OTFS_HPP_

#include <Module/Waveform/Ofdm.hpp>
#include <Tools/Math/matrix.h>
#include <fftw3.h>

namespace aff3ct
{
namespace module
{
    template <typename B = float>
    class Otfs: public Ofdm<B>
    {
    private:
        fftwf_complex *isfft_data_inplace;
        fftwf_plan isfft_doppler;
        fftwf_plan isfft_delay;

        fftwf_complex *sfft_data_inplace;
        fftwf_plan sfft_doppler;
        fftwf_plan sfft_delay;

    public:
        Otfs(int M, const std::vector<int> cp, int N);
        virtual ~Otfs();
        
        // data in X_K is column-major, which is M x N
        void modulate(const B *X_K, B *Y_K,int frame_id = -1); 
        void demodulate(const B *X_K, B *Y_K,int frame_id = -1); 
    protected:
        void _isfft(const B *X_K, B *Y_K); // IFFT for doppler(row), FFT for delay(col)
        void _sfft(const B *X_K, B *Y_K); // FFT for doppler(row), IFFT for delay(col)
    };
    
}
}
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Module/Waveform/OTFS/Otfs.hxx"
#endif

#endif