#ifndef FM_Index_misc_h
#define FM_Index_misc_h

template <class InputIterator, class Size, class OutputIterator, class UnaryPredicate>
OutputIterator copy_n_until(InputIterator first, Size n, OutputIterator result, UnaryPredicate pred)
{
    while (n>0) {
        if(pred(*first)) break;
        *result = *first;
        ++result; ++first;
        --n;
    }
    return result;
}

#endif
