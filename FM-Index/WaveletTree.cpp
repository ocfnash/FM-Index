#include <stdexcept>

#include "WaveletTree.h"
#include "serializing.h"

bool WaveletTree::is_leaf(void) const
{
    return left == nullptr; // Tree is full so could use right instead.
}

void WaveletTree::fill_alphabet(const char * s, const size_t len_s)
{
    std::set<unsigned char> a_set;
    for(size_t i = 0; i < len_s; i++)
    {
        a_set.insert(s[i]);
        if(a_set.size() == 256) break; // Minor optimisation.
    }
    alphabet_owner = std::unique_ptr<char[]>(new char[a_set.size()]);
    std::copy(a_set.cbegin(), a_set.cend(), alphabet_owner.get());
    alphabet_begin = alphabet_owner.get();
    alphabet_end = alphabet_begin + a_set.size();
}

bool WaveletTree::belongs_left(const char c) const
{
    return static_cast<unsigned char>(c) <=
           static_cast<unsigned char>(*(alphabet_begin + (1 + alphabet_end - alphabet_begin) / 2 - 1));
}

WaveletTree::WaveletTree(std::string & s,
                         const bool clear_s,
                         const char * alphabet_begin,
                         const char * alphabet_end)
    : alphabet_begin(alphabet_begin),
      alphabet_end(alphabet_end)
{
    if(s.size() == 0) throw std::length_error("Cannot construct zero-length WaveletTree");
    if(alphabet_begin == nullptr) fill_alphabet(s.c_str(), s.size());

    std::string s_left, s_right;
    std::vector<bool> data_v;
    for(std::string::iterator it = s.begin(); it != s.end(); it++)
    {
        if(belongs_left(*it))
        {
            s_left.push_back(*it);
            data_v.push_back(1);
        }
        else
        {
            s_right.push_back(*it);
            data_v.push_back(0);
        }
    }
    if(clear_s)
    {
        s.clear();
        s.shrink_to_fit();
    }
    data = std::unique_ptr<BitVector>(new BitVector(data_v));
    data_v.clear();
    data_v.shrink_to_fit();

    size_t alphabet_size = this->alphabet_end - this->alphabet_begin;
    const char *alphabet_mid = this->alphabet_begin + (1 + alphabet_size) / 2;
    if(alphabet_size <= 2)
    {
        left = nullptr;
        right = nullptr;
    }
    else
    {
        left = std::unique_ptr<WaveletTree>(new WaveletTree(s_left, true, this->alphabet_begin, alphabet_mid));
        right = std::unique_ptr<WaveletTree>(new WaveletTree(s_right, true, alphabet_mid, this->alphabet_end));
    }
}

size_t WaveletTree::size(void) const
{
    return data->size();
}

std::string WaveletTree::get_alphabet(void) const
{
    return std::string(alphabet_begin, alphabet_end);
}

size_t WaveletTree::cum_freq(const char c) const
{
    if(is_leaf())
    {
        if(belongs_left(c)) return 0;
        else return data->rank1(size() - 1);
    }
    else
    {
        if(belongs_left(c)) return left->cum_freq(c);
        else return left->size() + right->cum_freq(c);
    }
}

size_t WaveletTree::rank(const size_t i, const char c) const
{
    if(i >= data->size()) throw std::out_of_range("WaveletTree rank out of range");
    if(is_leaf())
    {
        if(c != *alphabet_begin && c != *(alphabet_end-1)) return 0; // c outside alphabet.
        if(belongs_left(c)) return data->rank1(i);
        else return data->rank0(i);
    }
    else
    {
        if(belongs_left(c)) return data->rank1(i) >= 1 ? left->rank(data->rank1(i)-1, c) : 0;
        else return data->rank0(i) >= 1 ? right->rank(data->rank0(i)-1, c) : 0;
    }
}

char WaveletTree::select(const size_t i) const
{
    if(i >= data->size()) throw std::out_of_range("WaveletTree select out of range");
    if(is_leaf()) return alphabet_begin[1-data->select(i)]; // NB: true == 1 is in C++ spec.
    else return data->select(i) == 1 ? left->select(data->rank1(i)-1) : right->select(data->rank0(i)-1);
}

WaveletTree::WaveletTree(std::istreambuf_iterator<char> serial_data)
{
    size_t alphabet_size;
    deserialize_from_chars(serial_data, alphabet_size);
    alphabet_owner = std::unique_ptr<char[]>(new char[alphabet_size]);
    for(size_t i = 0; i < alphabet_size; i++)
    {
        alphabet_owner[i] = *serial_data;
        ++serial_data;
    }
    alphabet_begin = alphabet_owner.get();
    alphabet_end = alphabet_begin + alphabet_size;
    data = std::unique_ptr<BitVector>(new BitVector(serial_data));
    char n_children = *serial_data;
    ++serial_data;
    if(n_children > 0)
    {
        //assert(*serial_data == 2);
        left = std::unique_ptr<WaveletTree>(new WaveletTree(serial_data));
        right = std::unique_ptr<WaveletTree>(new WaveletTree(serial_data));
    }
}

void WaveletTree::serialize(std::ostreambuf_iterator<char> serial_data) const
{
    /* Note that because we write out a copy of the alphabet, each child
       has its own copy of the alphabet when serialized. This is different
       to when a WaveletTree is constructed from a std::string when children
       will just have pointers to the parent alphabet (this saves children
       repeating work of parent scanning for alphabet). For this reason
       serializing and deserializing will not produce identical objects,
       though the difference is invisible to the user. */
    serialize_as_chars(serial_data, alphabet_end - alphabet_begin);
    for(char * p = const_cast<char *>(alphabet_begin); p != alphabet_end; p++)
    {
        *serial_data = *p;
        ++serial_data;
    }
    data->serialize(serial_data);
    if(left != nullptr)
    {
        //assert(right != nullptr);
        *serial_data = 2; // 2 children.
        ++serial_data;
        left->serialize(serial_data);
        right->serialize(serial_data);
    }
    else
    {
        *serial_data = 0; // 0 children.
        ++serial_data;
    }
}
