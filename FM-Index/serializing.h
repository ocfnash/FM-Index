#ifndef FM_Index_serializing_h
#define FM_Index_serializing_h

#include <sstream>
#include <iterator>

// I suspect these functions are VERY slow. I'll wait and see if
// speed becomes annoying before thinking about optimizing.

template <typename T>
void serialize_as_chars(std::ostreambuf_iterator<char> s, const T & x)
{
    const char * p = reinterpret_cast<const char *>(&x);
    for(size_t i = 0; i < sizeof(T); i++)
    {
        *s = *p;
        ++s; ++p;
    }
}

template <typename T>
void deserialize_from_chars(std::istreambuf_iterator<char> s, T & x)
{
    char * p = reinterpret_cast<char *>(&x);
    for(size_t i = 0; i < sizeof(T); i++)
    {
        *p = *s;
        ++s; ++p;
    }
}

#endif
