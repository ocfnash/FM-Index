#ifndef __FM_Index__FMIndex__
#define __FM_Index__FMIndex__

#include <map>
#include <list>
#include <iterator>

#include "WaveletTree.h"

class FMIndex
{
public:
    class const_iterator : public std::iterator<std::forward_iterator_tag, char>
    {
    private:
        const std::unique_ptr<WaveletTree> & BWT_or_BWTr;
        const size_t end_idx;
        const std::map<char, size_t> & C;
        size_t i; // A row index in the hypothetical matrix whose last column is BWT_or_BWTr.
        char c; // The character at the end of the row in the hypothetical matrix (i.e., in BWT_or_BWTr).

    public:
        const_iterator(const std::unique_ptr<WaveletTree> & BWT_or_BWTr,
                       const size_t end_idx,
                       const std::map<char, size_t> & C,
                       const size_t i);

        const_iterator(const const_iterator & it)
            : BWT_or_BWTr(it.BWT_or_BWTr),
              end_idx(it.end_idx),
              C(it.C),
              i(it.i),
              c(it.c) { }

        //const_iterator & operator=(const const_iterator & it);

        bool operator==(const const_iterator & it);

        bool operator!=(const const_iterator & it);

        const_iterator & operator++(void); // Prefix

        const_iterator operator++(int); // Postfix

        char operator*(void) const;

        //char * operator->(void) const;

        bool at_end(void) const;
    };

    typedef const_iterator const_reverse_iterator; // What is type-safe way of doing this?

private:
    std::unique_ptr<WaveletTree> BWT_as_wt, BWTr_as_wt;
    size_t BWT_end_idx, BWTr_end_idx;
    std::map<char, size_t> C;

    static size_t BWT_idx_from_row_idx(const size_t i, const size_t end_idx);

    template <typename ForwardIterator>
    std::pair<size_t, size_t> backward_search(ForwardIterator i_pattern,
                                              ForwardIterator i_pattern_end,
                                              const std::unique_ptr<WaveletTree> & BWT_or_BWTr,
                                              const size_t end_idx) const;

    void msd_sort(size_t * fr_perm,
                  const size_t len,
                  const_iterator ** text_iters,
                  const size_t depth_left) const;

    void populate_C(void);

public:
    FMIndex(const std::string & s);

    FMIndex(std::istreambuf_iterator<char> serial_data);

    size_t findn(const std::string & pattern) const;

    size_t find(std::list<std::pair<const_iterator, const_reverse_iterator>> & matches,
                const std::string & pattern,
                const size_t max_context = 100) const;

    std::list<std::string> find_lines(const std::string & pattern,
                                      const char new_line_char = '\n',
                                      const size_t max_context = 100) const;

    size_t size(void) const;

    const_iterator begin(void) const;

    void serialize(std::ostreambuf_iterator<char> serial_data) const;

    void serialize_to_file(const std::string & filename) const;

    static FMIndex * new_from_serialized_file(const std::string & filename);
};

#endif /* defined(__FM_Index__FMIndex__) */
