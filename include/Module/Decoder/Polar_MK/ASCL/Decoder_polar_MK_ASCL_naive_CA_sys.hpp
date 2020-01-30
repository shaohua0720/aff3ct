/*!
 * \file
 * \brief Class module::Decoder_polar_MK_ASCL_naive_CA_sys.
 */
#ifndef DECODER_POLAR_MK_ASCL_NAIVE_CA_SYS
#define DECODER_POLAR_MK_ASCL_NAIVE_CA_SYS

#include <functional>
#include <vector>
#include <memory>

#include "Tools/Code/Polar/Polar_code.hpp"
#include "Module/CRC/CRC.hpp"
#include "Module/Decoder/Polar_MK/SC/Decoder_polar_MK_SC_naive_sys.hpp"
#include "Module/Decoder/Polar_MK/SCL/CRC/Decoder_polar_MK_SCL_naive_CA_sys.hpp"

namespace aff3ct
{
namespace module
{
template <typename B = int, typename R = float>
class Decoder_polar_MK_ASCL_naive_CA_sys : public Decoder_polar_MK_SCL_naive_CA_sys<B,R>
{
private:
	std::shared_ptr<Decoder_polar_MK_SC_naive_sys<B,R>> sc_decoder;

public:
	Decoder_polar_MK_ASCL_naive_CA_sys(const int& K,
	                                   const int& N,
	                                   const int& L,
	                                   const tools::Polar_code& code,
	                                   const std::vector<bool>& frozen_bits,
	                                   const std::vector<std::vector<std::function<R(const std::vector<R> &LLRs,
	                                                                                 const std::vector<B> &bits)>>> &lambdas,
	                                   const CRC<B>& crc,
	                                   const int n_frames = 1);

	Decoder_polar_MK_ASCL_naive_CA_sys(const int& K,
	                                   const int& N,
	                                   const int& L,
	                                   const tools::Polar_code& code,
	                                   const std::vector<bool>& frozen_bits,
	                                   const CRC<B>& crc,
	                                   const int n_frames = 1);

	virtual ~Decoder_polar_MK_ASCL_naive_CA_sys() = default;

	virtual Decoder_polar_MK_ASCL_naive_CA_sys<B,R>* clone() const;

	virtual void notify_noise_update();

protected:
	virtual void deep_copy(const Decoder_polar_MK_ASCL_naive_CA_sys<B,R> &m);

	void _decode_siho   (const R *Y_N, B *V_K, const int frame_id);
	void _decode_siho_cw(const R *Y_N, B *V_N, const int frame_id);
};
}
}

#endif /* DECODER_POLAR_MK_ASCL_NAIVE_CA_SYS */