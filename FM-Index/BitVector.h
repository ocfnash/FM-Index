#ifndef __FM_Index__BitVector__
#define __FM_Index__BitVector__

#include <vector>

class BitVector
{
private:
    typedef uint64_t block_t; // NB: Type must match __builtin_popcount or __builtin_popcountl etc in rank method.
    static const size_t superblock_sz_bits = 512; // Must be multiple of size_of_data_t_bits.
    static const size_t size_of_data_t_bits = 8 * sizeof(block_t);
    std::unique_ptr<block_t[]> data; // Should be slightly faster than std::bitset<...> *data. (Also not sure size overhead of std::bitset<...> is 0.)
    std::unique_ptr<uint32_t[]> superblock_ranks; // 32 bits might not be enough but unlikely on "real" data. Overflow throws exception anyway.
    size_t q, r;

    size_t rank(const size_t i) const;

public:
    BitVector(const std::vector<bool> & data);

    BitVector(std::istreambuf_iterator<char> serial_data);

    size_t rank0(const size_t i) const;

    size_t rank1(const size_t i) const;

    bool select(const size_t i) const;

    size_t size(void) const;

    void serialize(std::ostreambuf_iterator<char> serial_data) const;
};

#endif /* defined(__FM_Index__BitVector__) */