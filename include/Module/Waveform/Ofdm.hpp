#ifndef OFDM_HPP_
#define OFDM_HPP_

#include <iostream>
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
                generate,
                SIZE
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
            Ofdm(const int M, const int N = 14);
            virtual ~Ofdm() = default;
            virtual Ofdm<B> *clone() const;

            template <class AB = std::allocator<B>>
            void modulate(const std::vector<B, AB> &X_K, std::vector<B, AB> &Y_K, int frame_id = -1);
            template <typename AB = std::allocator<B>>
            void demodulate(const std::vector<B, AB> &X_K, std::vector<B, AB> &Y_K, int frame_id = -1);
        };
    }

} // namespace aff3ct

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Module/Waveform/Ofdm.hxx"
#endif

#endif