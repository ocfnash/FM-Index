#include <algorithm>
#include <fstream>
#include <iterator>

#include "FMIndex.h"
#include "openbwt.h"
#include "serializing.h"
#include "misc.h"

FMIndex::const_iterator::const_iterator(const std::unique_ptr<WaveletTree> & BWT_or_BWTr,
                                        const size_t end_idx,
                                        const std::map<char, size_t> & C,
                                        const size_t i)
    : BWT_or_BWTr(BWT_or_BWTr),
      end_idx(end_idx),
      C(C),
      i(i)
{
    if(i > BWT_or_BWTr->size()) throw std::out_of_range("Attempt to create FMIndex::const_iterator with out-of-bounds index");
    if(!at_end()) c = BWT_or_BWTr->select(BWT_idx_from_row_idx(i, end_idx));
}

bool FMIndex::const_iterator::operator==(const FMIndex::const_iterator::const_iterator & it)
{
    return (BWT_or_BWTr == it.BWT_or_BWTr &&
            end_idx == it.end_idx &&
            C == it.C &&
            i == it.i);
}

bool FMIndex::const_iterator::operator!=(const FMIndex::const_iterator & it)
{
    return !operator==(it);
}

FMIndex::const_iterator & FMIndex::const_iterator::operator++(void)
{
    if(at_end()) throw std::overflow_error("Attempt to increment ended const_iterator");
    i = C.find(c)->second + BWT_or_BWTr->rank(BWT_idx_from_row_idx(i, end_idx), c);
    if(!at_end()) c = BWT_or_BWTr->select(BWT_idx_from_row_idx(i, end_idx));
    return *this;
}

FMIndex::const_iterator FMIndex::const_iterator::operator++(int)
{
    FMIndex::const_iterator it(*this);
    operator++();
    return it;
}

char FMIndex::const_iterator::operator*(void) const
{
    if(at_end()) throw std::overflow_error("Attempt to read character from ended const_iterator");
    return c;
}

bool FMIndex::const_iterator::at_end(void) const
{
    return i == end_idx;
}

size_t FMIndex::BWT_idx_from_row_idx(const size_t i, const size_t end_idx)
{
    if(i == end_idx) throw std::invalid_argument("No BWT index corresponding to row end_idx");
    return (i > end_idx ? i-1 : i);
}

template <typename ForwardIterator>
std::pair<size_t, size_t> FMIndex::backward_search(ForwardIterator i_pattern,
                                                   ForwardIterator i_pattern_end,
                                                   const std::unique_ptr<WaveletTree> & BWT_or_BWTr,
                                                   const size_t end_idx) const
{
    std::pair<size_t, size_t> no_matches(1, 0);

    std::map<char, size_t>::const_iterator C_it = C.find(*i_pattern);
    if(C_it == C.end()) return no_matches;
    /* lb and ub define the half-open interval (lb, ub] of row indexes into the
       hypothetical matrix for which the rows are prefixed by the pattern. */
    size_t lb = C_it->second;
    size_t ub = (++C_it == C.end() ? BWT_or_BWTr->size() : C_it->second);
    for(++i_pattern; i_pattern != i_pattern_end; ++i_pattern)
    {
        C_it = C.find(*i_pattern);
        if(C_it == C.end()) return no_matches;
        lb = C_it->second + BWT_or_BWTr->rank(BWT_idx_from_row_idx(lb, end_idx), *i_pattern);
        ub = C_it->second + BWT_or_BWTr->rank(BWT_idx_from_row_idx(ub, end_idx), *i_pattern);
        if(ub <= lb) return no_matches;
    }
    return std::make_tuple(lb + 1, ub + 1); // Return as more-conventional half-open interval [lb, ub)
}

void FMIndex::msd_sort(size_t * fr_perm,
                       const size_t len,
                       const_iterator ** text_iters,
                       const size_t depth_left) const
{
    // Hybrid of most-significant-digit radix sort and std::sort.
    // TODO: Replace std::sort with American Flag Sort.
    std::sort(fr_perm, fr_perm + len,
              [text_iters](const size_t i, const size_t j) -> bool
              {
                  if(text_iters[j]->at_end()) return false;
                  if(text_iters[i]->at_end()) return true;
                  return static_cast<unsigned char>(**text_iters[i]) <
                         static_cast<unsigned char>(**text_iters[j]);
              });
    int c = -1;
    if(!text_iters[fr_perm[0]]->at_end())
    {
        c = static_cast<unsigned char>(**text_iters[fr_perm[0]]);
        ++(*(text_iters[fr_perm[0]]));
    }
    size_t j = 0;
    for(size_t i = 1; i < len; i++)
    {
        int cc = -1;
        if(!text_iters[fr_perm[i]]->at_end())
        {
            cc = static_cast<unsigned char>(**text_iters[fr_perm[i]]);
            ++(*(text_iters[fr_perm[i]]));
        }
        if(cc != c)
        {
            if(i - j > 1 && depth_left > 0) msd_sort(fr_perm + j, i - j, text_iters, depth_left - 1);
            j = i;
            c = cc;
        }
    }
    if(len - j > 1 && depth_left > 0) msd_sort(fr_perm + j, len - j, text_iters, depth_left - 1);
}

void FMIndex::populate_C(void)
{
    std::string alphabet = BWT_as_wt->get_alphabet();
    for(std::string::iterator c = alphabet.begin(); c != alphabet.end(); c++)
        C[*c] = BWT_as_wt->cum_freq(*c);
}

FMIndex::FMIndex(const std::string & s)
{
    if(s.empty()) throw std::length_error("Cannot construct zero-length FMIndex");

    // Build BWT_as_wt:
    std::string s_BWT;
    s_BWT.resize(s.size());
    BWT_end_idx = BWT((const unsigned char *) s.c_str(), (unsigned char *) &s_BWT[0], (int) s.size()); // C-style casts for C-function BWT.
    BWT_as_wt = std::unique_ptr<WaveletTree>(new WaveletTree(s_BWT, false));

    // Build BWTr_as_wt:
    std::string s_rev(s.rbegin(), s.rend()); // Uugh, waste of time+space. Ideally should teach BWT to optionally sort right-left.
    BWTr_end_idx = BWT((const unsigned char *) s_rev.c_str(), (unsigned char *) &s_BWT[0], (int) s.size()); // C-style casts for C-function BWT.
    BWTr_as_wt = std::unique_ptr<WaveletTree>(new WaveletTree(s_BWT));

    populate_C();
}

size_t FMIndex::findn(const std::string & pattern) const
{
    if(pattern.empty()) throw std::length_error("Cannot search for zero-length pattern");

    size_t lb, ub;
    std::tie(lb, ub) = backward_search(pattern.rbegin(), pattern.rend(), BWT_as_wt, BWT_end_idx);
    return ub <= lb ? 0 : ub - lb;
}

size_t FMIndex::find(std::list<std::pair<const_iterator, const_reverse_iterator>> & matches,
                     const std::string & pattern,
                     const size_t max_context) const
{
    if(pattern.empty()) throw std::length_error("Cannot search for zero-length pattern");

    size_t lb, ub, lbr, ubr;
    std::tie(lb, ub) = backward_search(pattern.rbegin(), pattern.rend(), BWT_as_wt, BWT_end_idx);
    if(ub <= lb) return 0;
    std::tie(lbr, ubr) = backward_search(pattern.begin(), pattern.end(), BWTr_as_wt, BWTr_end_idx);
    assert(ub-lb == ubr-lbr);

    // Slightly painful last step: find permutation that matches up indices into BWT_as_wt and BWTr_as_wt.
    size_t n_matches = ub - lb;
    size_t * fr_perm = new size_t[n_matches];
    const_iterator ** text_iters;
    text_iters = new const_iterator*[n_matches];
    for(size_t i = 0; i < n_matches; i++)
    {
        text_iters[i] = new const_iterator(BWTr_as_wt, BWTr_end_idx, C, lbr + i);
        fr_perm[i] = i;
    }
    msd_sort(fr_perm, n_matches, text_iters, max_context);
    for(size_t i = 0; i < n_matches; i++)
    {
        matches.push_back(std::make_tuple(const_iterator(BWTr_as_wt, BWTr_end_idx, C, lbr + fr_perm[i]),
                                          const_reverse_iterator(BWT_as_wt, BWT_end_idx, C, lb + i)));
        delete text_iters[i];
    }
    delete [] text_iters;
    delete [] fr_perm;

    return n_matches;
}

std::list<std::string> FMIndex::find_lines(const std::string & pattern,
                                           const char new_line_char,
                                           const size_t max_context) const
{
    std::list<std::pair<const_iterator, const_reverse_iterator>> matches;
    find(matches, pattern, max_context);

    std::list<std::string> l;
    for(auto & match : matches)
    {
        std::ostringstream context_before_ss;
        try
        {
            copy_n_until(match.second,
                         max_context,
                         std::ostream_iterator<char>(context_before_ss),
                         [new_line_char](char c) -> bool
                         {
                             return c == new_line_char;
                         });
        }
        catch(std::overflow_error & e) { };
        std::ostringstream context_after_ss;
        try
        {
            copy_n_until(match.first,
                         max_context,
                         std::ostream_iterator<char>(context_after_ss),
                         [new_line_char](char c) -> bool
                         {
                             return c == new_line_char;
                         });
        }
        catch(std::overflow_error & e) { };
        std::string context_before_s = context_before_ss.str();
        l.push_back(std::string(context_before_s.rbegin(), context_before_s.rend()) +
                    pattern +
                    context_after_ss.str());
    }

    return l;
}

size_t FMIndex::size(void) const
{
    //assert(BWT_as_wt->size() == BWTr_as_wt->size());
    return BWT_as_wt->size();
}

FMIndex::const_iterator FMIndex::begin(void) const
{
    return const_iterator(BWTr_as_wt, BWTr_end_idx, C, 0);
}

void FMIndex::serialize_to_file(const std::string & filename) const
{
    /* Intended for use in python (via Cython wrapper). Important as
       will save copying *large* amounts of data from C++ to python */
    std::ofstream f{filename, std::ios::binary};
    serialize(std::ostreambuf_iterator<char>{f});
}

FMIndex * FMIndex::new_from_serialized_file(const std::string & filename)
{
    std::ifstream f{filename, std::ios::binary};
    return new FMIndex{std::istreambuf_iterator<char>{f}};
}

FMIndex::FMIndex(std::istreambuf_iterator<char> serial_data)
{
    BWT_as_wt = std::unique_ptr<WaveletTree>(new WaveletTree(serial_data));
    deserialize_from_chars(serial_data, BWT_end_idx);
    BWTr_as_wt = std::unique_ptr<WaveletTree>(new WaveletTree(serial_data));
    deserialize_from_chars(serial_data, BWTr_end_idx);
    populate_C();
}

void FMIndex::serialize(std::ostreambuf_iterator<char> serial_data) const
{
    BWT_as_wt->serialize(serial_data);
    serialize_as_chars(serial_data, BWT_end_idx);
    BWTr_as_wt->serialize(serial_data);
    serialize_as_chars(serial_data, BWTr_end_idx);
}
