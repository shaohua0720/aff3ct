#include "Tools/Math/lconv.h"
#include <vector>

namespace aff3ct
{
    namespace tools
    {
        template <typename T, class AT = std::allocator<T>>
        void conv(const std::vector<T, AT> seq1, const int len1,
                  const std::vector<T, AT> seq2, const int len2,
                  std::vector<T, AT> &out, int &out_len)
        {
            out_len = len1 + len2 - 1;
            if (out.size()<out_len)
            {
                out.resize(out_len);
            }
            
            for (int i = 0; i < len1; i++)
                for (int j = 0; j < len2; j++)
                {
                    out[i + j] = out[i + j] + seq1[i] * seq2[j];
                }
        }

    }
}