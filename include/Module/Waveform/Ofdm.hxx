#include "Tools/Exception/exception.hpp"
#include "Module/Waveform/Ofdm.hpp"
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
        Ofdm<B>::Ofdm(const int M, const int N) : Module(), M(M), N(N)
        {
            const std::string name = "Waveform";
            this->set_name(name);
            this->set_short_name(name);

            // create related modulate and demodulates
            auto& p = this->create_task("modulate");
            auto ps_X_K = this->template create_socket_in<B>(p, "X_K", 1);
            this->create_codelet(p, [ps_X_K](Module &m, Task &t, const size_t frame_id) -> int
            {
                auto &ofdm = static_cast<Ofdm<B>&>(m);
            }
            

        }

        template <typename B>
        Ofdm<B> *Ofdm<B>::clone() const
        {
            throw tools::unimplemented_error(__FILE__, __LINE__, __func__);
        }
    }
}