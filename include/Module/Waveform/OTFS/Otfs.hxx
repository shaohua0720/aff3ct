#include "Module/Waveform/OTFS/Otfs.hpp"
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
        }
        
        template <typename B>
        Otfs<B>::~Otfs()
        {

        }
        template <typename B>
        void Otfs<B>::modulate(const B *X_K, B *Y_K,int frame_id)
        {

        }
        template <typename B>
        void Otfs<B>::demodulate(const B *X_K, B *Y_K,int frame_id)
        {

        }
    }

}