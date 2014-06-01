#ifndef __FM_Index__WaveletTree__
#define __FM_Index__WaveletTree__

#include <set>
#include <string>

#include "BitVector.h"

class WaveletTree
{
private:
    std::unique_ptr<char[]> alphabet_owner;
    const char *alphabet_begin, *alphabet_end;
    std::unique_ptr<BitVector> data;
    std::unique_ptr<WaveletTree> left, right;

    bool is_leaf(void) const;

    bool belongs_left(const char c) const;

    void fill_alphabet(const char * s, const size_t len_s);

public:
    WaveletTree(std::string & s,
                const bool clear_s = true,
                const char * alphabet_begin = nullptr,
                const char * alphabet_end = nullptr);

    WaveletTree(std::istreambuf_iterator<char> serial_data);

    size_t size(void) const;

    std::string get_alphabet(void) const;

    size_t cum_freq(const char c) const;

    size_t rank(const size_t i, const char c) const;

    char select(const size_t i) const;

    void serialize(std::ostreambuf_iterator<char> serial_data) const;
};

#endif /* defined(__FM_Index__WaveletTree__) */
