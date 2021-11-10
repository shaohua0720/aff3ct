#include <exception>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <mipp.h>
#ifdef AFF3CT_MPI
#include <mpi.h>
#endif
#include <cli.hpp>

#include "Tools/types.h"
#include "Tools/version.h"
#include "Tools/Display/rang_format/rang_format.h"
#include "Launcher/Launcher.hpp"
#include "Factory/Launcher/Launcher.hpp"
#include "aff3ct.hpp"
#include <complex>

#include "aff3ct.hpp"

using namespace aff3ct;
using namespace module;

void print_version()
{
#if defined(_WIN64) || defined(_WIN32) || defined(WIN32) || defined(__MINGW32__)
	std::string os = "Windows";
#elif defined(__linux__) || defined(__linux)
	std::string os = "Linux";
#elif defined(__APPLE__)
	std::string os = "Mac OS X";
#else
	std::string os = "Unknown OS";
#endif

#if defined(__x86_64) || defined(__x86_64__) || defined(_WIN64) || defined(__aarch64__)
	std::string prec = " 64-bit";
#elif defined(__i386) || defined(__i386__) || defined(i386) || defined(_WIN32) || defined(__ARM_ARCH_7__)
	std::string prec = " 32-bit";
#else
	std::string prec = "";
#endif

#if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC)
	std::string compiler = "icpc";
#if defined(__INTEL_COMPILER)
	std::string compiler_version = std::to_string(__INTEL_COMPILER);
#elif defined(__ICL)
	std::string compiler_version = std::to_string(__ICL);
#else
	std::string compiler_version = std::to_string(__ICC);
#endif
	compiler_version = compiler_version.substr(0, 2) + "." + compiler_version.substr(2, compiler_version.size());
#elif defined(__clang__) || defined(__llvm__)
	std::string compiler = "clang++";
	std::string compiler_version = std::to_string(__clang_major__);
	compiler_version += "." + std::to_string(__clang_minor__);
#elif defined(__GNUG__) || (defined(__GNUC__) && defined(__cplusplus))
	std::string compiler = "g++";
	std::string compiler_version = std::to_string(__GNUC__);
	compiler_version += "." + std::to_string(__GNUC_MINOR__);
#elif defined(_MSC_VER)
	std::string compiler = "MSVC++";
	std::string compiler_version = std::to_string(_MSC_VER);
#else
	std::string compiler = "Unknown compiler";
	std::string compiler_version = "";
#endif
	std::string affect_version = tools::version() == "GIT-NOTFOUND" ? "" : tools::version();

#if defined(AFF3CT_MULTI_PREC)
	std::string precision = "8/16/32/64-bit";
#elif defined(AFF3CT_8BIT_PREC)
	std::string precision = "8-bit";
#elif defined(AFF3CT_16BIT_PREC)
	std::string precision = "16-bit";
#elif defined(AFF3CT_32BIT_PREC)
	std::string precision = "32-bit";
#elif defined(AFF3CT_64BIT_PREC)
	std::string precision = "64-bit";
#else
	std::string precision = "Unknown";
#endif

#if defined(AFF3CT_POLAR_BIT_PACKING)
	std::string bit_packing = "on";
#else
	std::string bit_packing = "off";
#endif

#if defined(AFF3CT_POLAR_BOUNDS)
	std::string polar_bounds = "on";
#else
	std::string polar_bounds = "off";
#endif

#if defined(AFF3CT_COLORS)
	std::string terminal_colors = "on";
#else
	std::string terminal_colors = "off";
#endif

#if defined(AFF3CT_BACKTRACE)
	std::string backtrace = "on";
#else
	std::string backtrace = "off";
#endif

#if defined(AFF3CT_EXT_STRINGS)
	std::string ext_strings = "on";
#else
	std::string ext_strings = "off";
#endif

#if defined(AFF3CT_MPI)
	std::string mpi = "on";
#else
	std::string mpi = "off";
#endif

#if defined(AFF3CT_CHANNEL_GSL)
	std::string gsl = "on";
#else
	std::string gsl = "off";
#endif

#if defined(AFF3CT_CHANNEL_MKL)
	std::string mkl = "on";
#else
	std::string mkl = "off";
#endif

	std::cout << "aff3ct (" << os << prec << ", " << compiler << "-" << compiler_version << ", "
			  << mipp::InstructionFullType << ") " << affect_version << std::endl;
	std::cout << "Compilation options:" << std::endl;
	std::cout << "  - Precision:         " << precision << std::endl;
	std::cout << "  - Polar bit packing: " << bit_packing << std::endl;
	std::cout << "  - Polar bounds:      " << polar_bounds << std::endl;
	std::cout << "  - Terminal colors:   " << terminal_colors << std::endl;
	std::cout << "  - Backtrace:         " << backtrace << std::endl;
	std::cout << "  - External strings:  " << ext_strings << std::endl;
	std::cout << "  - MPI:               " << mpi << std::endl;
	std::cout << "  - GSL:               " << gsl << std::endl;
	std::cout << "  - MKL:               " << mkl << std::endl;
	std::cout << "Copyright (c) 2016-2021 - MIT license." << std::endl;
	std::cout << "This is free software; see the source for copying conditions.  There is NO" << std::endl;
	std::cout << "warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << std::endl;
	exit(EXIT_SUCCESS);
}

int read_arguments(const int argc, const char **argv, factory::Launcher &params)
{
	cli::Argument_handler ah(argc, argv);

	cli::Argument_map_info args;
	cli::Argument_map_group arg_group;

	std::vector<std::string> cmd_warn, cmd_error;

	params.get_description(args);

	auto arg_vals = ah.parse_arguments(args, cmd_warn, cmd_error);

	bool display_help = false;
	try
	{
		params.store(arg_vals);
		ah.set_help_display_keys(params.display_keys);
	}
	catch (std::exception &)
	{
		display_help = true;
	}

	if (params.display_version)
		print_version();

	if (cmd_error.size() || display_help)
	{
		arg_group["sim"] = "Simulation parameter(s)";
		ah.print_help(args, arg_group, params.display_adv_help);

		if (cmd_error.size())
			std::cerr << std::endl;
		for (auto w = 0; w < (int)cmd_error.size(); w++)
			std::cerr << rang::tag::error << cmd_error[w] << std::endl;

		if (cmd_warn.size())
			std::cerr << std::endl;
		for (auto w = 0; w < (int)cmd_warn.size(); w++)
			std::cerr << rang::tag::warning << cmd_warn[w] << std::endl;
	}
	return (cmd_error.size() || display_help) ? EXIT_FAILURE : EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
	// std::vector<int> a = {1, 2, 3, 4, 5};
	// std::vector<int> b = {4, 5, 6};
	// std::vector<int> c;
	// int len = 0;
	// //c.resize(8);
	// tools::conv<int>(a, 5, b, 3, c, len);
	// std::cout << "linear convolution:";
	// for (auto ele : c)
	// {
	// 	std::cout << ele << " ";
	// }
	// std::cout << std::endl;
	// std::cout << "len:" << len << std::endl;

	// return 0;

	/******************************above code for test*************************************/

	// std::vector<int> cp = {1,2,1,1};
	// Ofdm<double> * ofdm = new Ofdm<double>(4,cp,4);

	// std::cout<<"SRC:"<<std::endl;
	// double a[32] = {
	// 	1.8033 , 1.7566,  0.1786 , 0.1893  , 2.9690 , 1.6081 ,  1.9659 , 1.6188,
	// 	0.2790 , 0.4408 ,  2.7060 , 2.4935  , 2.3315 , 3.8310 ,  0.9028 , 3.5956,
	// 	2.2513 , 2.1817 ,  0.2105 , 3.8126  , 2.3930 , 0.3826 ,  0.8399 , 4.3001,
	// 	0.5263 , 2.2085 ,  1.9261 , 4.9663  , 2.5478 , 0.2753 ,  1.5197 , 3.6756
	// 	};

	// for (int i = 0; i < 32; )
	// {
	// 	std::cout<<"("<<a[i]<<","<<a[i+1]<<")"<<" ";
	// 	i+=2;
	// }
	// std::cout<<std::endl;

	// // std::vector<int> cp = {0,0,0,0,};
	// // ofdm->setCPlength(cp);

	// double x_out[42];
	// ofdm->modulate(a,x_out,1);

	// std::cout<<"IFFT:"<<std::endl;
	// for(size_t i = 0; i < 42;)
	// {
	// 	std::cout<<"("<<x_out[i]/2<<","<<x_out[i+1]/2<<")"<<" ";
	// 	i+=2;
	// }
	// std::cout<<std::endl;
	// // double b[32] = {
	// // 	6.9168,5.1728,-2.5952,1.9358, 2.6278,1.5566,0.2638,-1.6388,
	// // 	6.2193,10.3609,-3.1546,-5.1934,-0.9983,-1.8173,-0.9504,-1.587,
	// // 	5.6947,10.677,-0.6292,2.4285,3.5939,-5.5484,0.3458,1.1697,
	// // 	6.5199,11.1257,-0.7308,1.5268,-0.3717,-6.1581,-3.3122,2.3396
	// // };
	// double c_out[32];
	// ofdm->demodulate(x_out,c_out,1);

	// std::cout<<"FFT:"<<std::endl;
	// for(size_t i = 0; i < 32; )
	// {
	// 	std::cout<<"("<<c_out[i]/4<<","<<c_out[i+1]/4<<")"<<" ";
	// 	i+=2;
	// }
	// std::cout<<std::endl;
	// return 0;
	/******************************** above for test **************************/
	// double a[32] = {
	// 	1.8033, 1.7566, 0.1786, 0.1893, 2.9690, 1.6081, 1.9659, 1.6188,
	// 	0.2790, 0.4408, 2.7060, 2.4935, 2.3315, 3.8310, 0.9028, 3.5956,
	// 	2.2513, 2.1817, 0.2105, 3.8126, 2.3930, 0.3826, 0.8399, 4.3001,
	// 	0.5263, 2.2085, 1.9261, 4.9663, 2.5478, 0.2753, 1.5197, 3.6756};
	float a[32] = {
		3.5468 , 4.4545,   3.2755 , 0.7465,   4.7987 , 4.0714,   3.7563 , 0.9830,
		3.7734 , 4.7965,   0.8131 , 1.2875,   1.7019 , 1.2176,   1.2755 , 1.2554,
		1.3801 , 2.7361,   0.5950 , 4.2036,   2.9263 , 4.6463,   2.5298 , 3.0802,
		3.3985 , 0.6931,   2.4918 , 1.2714,   1.1191 , 1.7499,   3.4954 , 2.3664
	};
	
	std::cout << "SRC:"<<std::endl;
	for (int i = 0; i < 32; )
	{
		std::cout<<"("<<a[i]<<","<<a[i+1]<<")"<<" ";
		i+=2;
	}

	std::vector<int> cp = {1, 2, 2, 1};
	Otfs<float> *otfs = new Otfs<float>(4, cp, 4);
	float b[44];
	otfs->modulate(a, b);

	std::cout<<std::endl;
	std::cout<<"Modulated:"<<std::endl;
	for(size_t i = 0; i < 44;)
	{
		std::cout<<"("<<b[i]<<","<<b[i+1]<<")"<<" ";
		i+=2;
	}
	float c[32];
	otfs->demodulate(b,c);

	std::cout<<std::endl;
	std::cout<<"Demodulated"<<std::endl;	
	for(size_t i = 0; i < 32;)
	{
		std::cout<<"("<<c[i]/64<<","<<c[i+1]/64<<")"<<" ";
		i+=2;
	}

	return 0;

	/******************************** above for test **************************/
	int exit_code = EXIT_SUCCESS;
#ifdef AFF3CT_MPI
	MPI_Init(nullptr, nullptr);
#endif

	factory::Launcher params("sim");
	if (read_arguments(argc, (const char **)argv, params) == EXIT_FAILURE)
		return EXIT_FAILURE;

	try
	{
		launcher::Launcher *launcher;
#ifdef AFF3CT_MULTI_PREC
		switch (params.sim_prec)
		{
		case 8:
			launcher = params.build<B_8, R_8, Q_8>(argc, (const char **)argv);
			break;
		case 16:
			launcher = params.build<B_16, R_16, Q_16>(argc, (const char **)argv);
			break;
		case 32:
			launcher = params.build<B_32, R_32, Q_32>(argc, (const char **)argv);
			break;
		case 64:
			launcher = params.build<B_64, R_64, Q_64>(argc, (const char **)argv);
			break;
		default:
			launcher = nullptr;
			break;
		}
#else
		launcher = params.build<B, R, Q>(argc, (const char **)argv);
#endif
		if (launcher != nullptr)
		{
			exit_code = launcher->launch();
			delete launcher;
		}
	}
	catch (std::exception const &e)
	{
		rang::format_on_each_line(std::cerr, std::string(e.what()) + "\n", rang::tag::error);
	}

#ifdef AFF3CT_MPI
	MPI_Finalize();
#endif

	return exit_code;
}