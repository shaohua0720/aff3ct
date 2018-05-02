#ifndef TERMINAL_EXIT_HPP_
#define TERMINAL_EXIT_HPP_

#include <string>
#include <chrono>
#include <iostream>

#include "Module/Monitor/EXIT/Monitor_EXIT.hpp"
#include "Tools/Noise/Noise.hpp"

#include "../Terminal.hpp"

namespace aff3ct
{
namespace tools
{
template <typename B = int, typename R = float>
class Terminal_EXIT : public Terminal
{
protected:
	const module::Monitor_EXIT<B,R> &monitor;
	std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds>  t_snr;
	unsigned short real_time_state;
	float   sig_a;
	Noise<float>* n;

public:
	explicit Terminal_EXIT(const module::Monitor_EXIT<B,R> &monitor);

	virtual ~Terminal_EXIT();

	void set_sig_a(float sig_a);

	void set_noise(const Noise<float>& noise);
	void set_noise(const Noise<double>& noise);

	void legend      (std::ostream &stream = std::cout);
	void temp_report (std::ostream &stream = std::cout);
	void final_report(std::ostream &stream = std::cout);

protected:
	void _report(std::ostream &stream);
};
}
}

#endif /* TERMINAL_EXIT_HPP_ */
