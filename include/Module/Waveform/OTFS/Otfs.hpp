#ifndef OTFS_HPP_
#define OTFS_HPP_

#include <Module/Waveform/Ofdm.hpp>

namespace aff3ct
{
namespace module
{
    template <typename B = double>
    class Otfs: public Ofdm<B>
    {
    private:
    public:
        Otfs(int M, const std::vector<int> cp, int N);
        virtual ~Otfs();
        void modulate(const B *X_K, B *Y_K,int frame_id = -1);
        void demodulate(const B *X_K, B *Y_K,int frame_id = -1);
    };
    
}
}
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Module/Waveform/OTFS/Otfs.hxx"
#endif

#endif