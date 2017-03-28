#include <stdexcept>

#include "matrix.h"

namespace aff3ct
{
namespace tools
{
template <typename T> 
inline void rgemm(const int M, const int N, const int K,
                  const mipp::vector<T> &A,
                  const mipp::vector<T> &tB,
                        mipp::vector<T> &tC)
{
	if (A.size() != unsigned(M * K))
		throw std::length_error("aff3ct::tools::rgemm: \"A.size()\" has to be equal to \"M\" * \"K\".");
	if (tB.size() != unsigned(K * N))
		throw std::length_error("aff3ct::tools::rgemm: \"tB.size()\" has to be equal to \"K\" * \"N\".");
	if (tC.size() != unsigned(M * N))
		throw std::length_error("aff3ct::tools::rgemm: \"tC.size()\" has to be equal to \"M\" * \"N\".");

	for (auto i = 0; i < M; i++)
		for (auto j = 0; j < N; j++)
		{
			T sum_r = 0;
			for (auto k = 0; k < K; k++)
				sum_r += A[i * K + k] * tB[j * K + k];

			tC[j * M + i] = sum_r;
		}
}

template <typename T>
inline void cgemm(const int M, const int N, const int K, 
                  const mipp::vector<T> &A, 
                  const mipp::vector<T> &tB, 
                        mipp::vector<T> &tC)
{
	if (A.size() != unsigned(M * K * 2))
		throw std::length_error("aff3ct::tools::cgemm: \"A.size()\" has to be equal to \"M\" * \"K\" * \"2\".");
	if (tB.size() != unsigned(K * N * 2))
		throw std::length_error("aff3ct::tools::cgemm: \"tB.size()\" has to be equal to \"K\" * \"N\" * \"2\".");
	if (tC.size() != unsigned(M * N * 2))
		throw std::length_error("aff3ct::tools::cgemm: \"tC.size()\" has to be equal to \"M\" * \"N\" * \"2\".");

	const T*  A_real =  A.data();
	const T*  A_imag =  A.data() + ( A.size() >> 1);
	const T* tB_real = tB.data();
	const T* tB_imag = tB.data() + (tB.size() >> 1);
	      T* tC_real = tC.data();
	      T* tC_imag = tC.data() + (tC.size() >> 1);

	for (auto i = 0; i < M; i++) 
	{
		for (auto j = 0; j < N; j++) 
		{
			T sum_r = 0, sum_i = 0;
			for (auto k = 0; k < K; k++) 
			{
				sum_r += A_real[i * K + k] * tB_real[j * K + k] - A_imag[i * K + k] * tB_imag[j * K + k];
				sum_i += A_imag[i * K + k] * tB_real[j * K + k] + A_real[i * K + k] * tB_imag[j * K + k];
			}

			tC_real[j * M + i] = sum_r;
			tC_imag[j * M + i] = sum_i;
		}
	}
}

template <typename T> 
inline void cgemm_r(const int M, const int N, const int K, 
                    const mipp::vector<T> &A, 
                    const mipp::vector<T> &tB, 
                          mipp::vector<T> &tC)
{
	if (A.size() != unsigned(M * K * 2))
		throw std::length_error("aff3ct::tools::cgemm_r: \"A.size()\" has to be equal to \"M\" * \"K\" * \"2\".");
	if (tB.size() != unsigned(K * N * 2))
		throw std::length_error("aff3ct::tools::cgemm_r: \"tB.size()\" has to be equal to \"K\" * \"N\" * \"2\".");
	if (tC.size() != unsigned(M * N * 1)) // because we only store the real part
		throw std::length_error("aff3ct::tools::cgemm_r: \"tC.size()\" has to be equal to \"M\" * \"N\" * \"1\".");

	const T*  A_real =  A.data();
	const T*  A_imag =  A.data() + ( A.size() >> 1);
	const T* tB_real = tB.data();
	const T* tB_imag = tB.data() + (tB.size() >> 1);
	      T* tC_real = tC.data();

	for (auto i = 0; i < M; i++) 
	{
		for (auto j = 0; j < N; j++) 
		{
			T sum_r = 0;
			for (auto k = 0; k < K; k++) 
				sum_r += A_real[i * K + k] * tB_real[j * K + k] - A_imag[i * K + k] * tB_imag[j * K + k];

			tC_real[j * M + i] = sum_r;
		}
	}
}

template <typename T>
inline void real_transpose(const int M, const int N,
                           const mipp::vector<T> &A,
                                 mipp::vector<T> &B)
{
	if (A.size() != unsigned(M * N))
		throw std::length_error("aff3ct::tools::real_transpose: \"A.size()\" has to be equal to \"M\" * \"N\".");
	if (B.size() != unsigned(N * M))
		throw std::length_error("aff3ct::tools::real_transpose: \"B.size()\" has to be equal to \"N\" * \"M\".");

	for (auto i = 0; i < M; i++)
		for (auto j = 0; j < N; j++)
			B[j*M+i] =  A[i*N+j];
}

template <typename T>
inline void complex_transpose(const int M, const int N,
                              const mipp::vector<T> &A,
                                    mipp::vector<T> &B)
{
	if (A.size() != unsigned(M * N * 2))
		throw std::length_error("aff3ct::tools::complex_transpose: \"A.size()\" has to be equal to "
		                        "\"M\" * \"N\" * \"2\".");
	if (B.size() != unsigned(N * M * 2))
		throw std::length_error("aff3ct::tools::complex_transpose: \"B.size()\" has to be equal to "
		                        "\"N\" * \"M\" * \"2\".");

	const T* A_real = A.data();
	const T* A_imag = A.data() + ( A.size() >> 1);
	      T* B_real = B.data();
	      T* B_imag = B.data() + (B.size() >> 1);

	for (auto i = 0; i < M; i++)
	{
		for (auto j = 0; j < N; j++)
		{
			B_real[j*M+i] =  A_real[i*N+j];
			B_imag[j*M+i] = -A_imag[i*N+j];
		}
	}
}
}
}
