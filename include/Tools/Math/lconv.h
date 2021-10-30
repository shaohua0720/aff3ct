#ifndef LCONV_H
#define LCONV_H
#include <vector>

namespace aff3ct
{
    namespace tools
    {
        template <typename T, class AT = std::allocator<T>>
        inline void conv(const std::vector<T, AT> seq1, const int len1,
                         const std::vector<T, AT> seq2, const int len2,
                         std::vector<T, AT> &out, int &out_len);
    }
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Math/lconv.hxx"
#endif
#endif