#include <stdexcept>

#include "BitVector.h"
#include "serializing.h"

BitVector::BitVector(const std::vector<bool> & data)
{
    if(data.size() == 0) throw std::length_error("Cannot construct zero-length BitVector");

    q = data.size() / superblock_sz_bits;
    r = data.size() % superblock_sz_bits;

    superblock_ranks = std::unique_ptr<uint32_t[]>(new uint32_t[q]);
    this->data = std::unique_ptr<block_t[]>(new block_t[1 + data.size() / size_of_data_t_bits]);

    size_t rk = 0;
    size_t i = 0;
    size_t j = 0;
    block_t x = 0;
    for(std::vector<bool>::const_iterator i_data = data.begin(); i_data != data.end(); i++, i_data++)
    {
        if(i % superblock_sz_bits == 0 && i > 0)
        {
            if(rk >= (1LU << 32)) throw std::overflow_error("Overflow error in superblock ranks.");
            superblock_ranks[i / superblock_sz_bits - 1] = static_cast<uint32_t>(rk);
        }
        if(i % size_of_data_t_bits == 0 && i > 0)
        {
            this->data[j++] = x;
            x = 0;
        }
        rk += *i_data;
        x = (x << 1) + *i_data;
    }
    if(i % size_of_data_t_bits != 0) x <<= size_of_data_t_bits - i % size_of_data_t_bits;
    this->data[j] = x;
}

size_t BitVector::rank(const size_t i) const
{
    if(i >= size()) throw std::out_of_range("BitVector rank out of range");

    /* We could get an on-average 2x speedup by scanning back from the
       next biggest superblock when we're nearer the end of the superblock
       than the beginning. Not bothering to implement for now as not crazy
       about speed issues at the moment. */
    size_t qq = i / superblock_sz_bits;
    size_t i1 = qq * (superblock_sz_bits / size_of_data_t_bits);
    size_t i2 = i / size_of_data_t_bits;
    size_t rr = i % size_of_data_t_bits;
    size_t rk = qq > 0 ? superblock_ranks[qq-1] : 0;
    // Slightly naughty using __builtin_popcountl below but let's live dangerously!
    for(size_t j = i1; j < i2; j++) rk += __builtin_popcountl(data[j]);
    return rk + __builtin_popcountl(data[i2] >> (size_of_data_t_bits - rr - 1));
}

size_t BitVector::rank0(const size_t i) const
{
    return i+1 - rank(i);
}

size_t BitVector::rank1(const size_t i) const
{
    return rank(i);
}

bool BitVector::select(const size_t i) const
{
    if(i >= size()) throw std::out_of_range("BitVector select out of range");

    size_t qq = i / size_of_data_t_bits;
    size_t rr = i % size_of_data_t_bits;
    return ((data[qq] >> (size_of_data_t_bits - rr - 1)) & 1) == 1;
}

size_t BitVector::size(void) const
{
    return r + q * superblock_sz_bits;
}

BitVector::BitVector(std::istreambuf_iterator<char> serial_data)
{
    size_t check;
    deserialize_from_chars(serial_data, check);
    if(check != superblock_sz_bits) throw std::runtime_error("Data for BitVector serialization has wrong superblock_sz_bits");
    deserialize_from_chars(serial_data, check);
    if(check != size_of_data_t_bits) throw std::runtime_error("Data for BitVector serialization has wrong size_of_data_t_bits");

    deserialize_from_chars(serial_data, q);
    deserialize_from_chars(serial_data, r);
    data = std::unique_ptr<block_t[]>(new block_t[1 + size() / size_of_data_t_bits]);
    for(size_t i = 0; i <= size() / size_of_data_t_bits; i++)
        deserialize_from_chars(serial_data, data[i]);
    superblock_ranks = std::unique_ptr<uint32_t[]>(new uint32_t[q]);
    for(size_t i = 0; i < q; i++)
        deserialize_from_chars(serial_data, superblock_ranks[i]);
}

void BitVector::serialize(std::ostreambuf_iterator<char> serial_data) const
{
    // I do not understand why the below casts are needed but I get a
    // *linker* error without them!
    serialize_as_chars(serial_data, reinterpret_cast<size_t>(superblock_sz_bits));
    serialize_as_chars(serial_data, reinterpret_cast<size_t>(size_of_data_t_bits));
    serialize_as_chars(serial_data, q);
    serialize_as_chars(serial_data, r);
    for(size_t i = 0; i <= size() / size_of_data_t_bits; i++)
        serialize_as_chars(serial_data, data[i]);
    for(size_t i = 0; i < q; i++)
        serialize_as_chars(serial_data, superblock_ranks[i]);
}
