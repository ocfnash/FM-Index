#include <string>
#include <sstream>

#include "gtest/gtest.h"
#include "BitVector.h"
#include "WaveletTree.h"
#include "FMIndex.h"
#include "openbwt.h"
#include "serializing.h"

class BitVectorTest : public ::testing::Test
{
protected:
    BitVector *zero_bv,
              *one_bv,
              *random1_bv,
              *random2_bv,
              *random3_bv;
    const std::vector<bool> random1_v{1,1,1,0,0,0,1,0,0,1},
                            random2_v{0,1,1,0,1,0,1,0},
                            random3_v{0,0,0,0,0,1,1,1,1,1,
                                      0,1,0,1,0,1,0,1,1,0,
                                      1,0,1,0,1,0,1,1,1,0,
                                      1,0,0,1,0,1,0,1,0,1,
                                      1,1,1,0,0,0,0,1,0,1,
                                      0,0,0,0,0,0,0,0,0,0,
                                      1,1,1,0,1,0,1,0,1,1};

    virtual void SetUp(void)
    {
        zero_bv = new BitVector(std::vector<bool>{0});
        one_bv = new BitVector(std::vector<bool>{1});
        random1_bv = new BitVector(random1_v);
        random2_bv = new BitVector(random2_v);
        random3_bv = new BitVector(random3_v);
    }

    virtual void TearDown(void)
    {
        delete zero_bv;
        delete one_bv;
        delete random1_bv;
        delete random2_bv;
        delete random3_bv;
    }
};

TEST(BWT, Mississippi)
{
    std::string s("mississippi");
    char t[s.size()];
    int pidx = BWT(reinterpret_cast<const unsigned char *>(s.c_str()),
                   reinterpret_cast<unsigned char *>(t),
                   static_cast<int>(s.size())); // Sorry about the casting.
    ASSERT_EQ(5, pidx);
    ASSERT_EQ(std::string("ipssmpissii"), std::string(t, s.size()));
}

TEST_F(BitVectorTest, Size)
{
    ASSERT_EQ(1, zero_bv->size());
    ASSERT_EQ(1, one_bv->size());
    ASSERT_EQ(random1_v.size(), random1_bv->size());
    ASSERT_EQ(random2_v.size(), random2_bv->size());
    ASSERT_EQ(random3_v.size(), random3_bv->size());
}

TEST_F(BitVectorTest, Rank)
{
    ASSERT_EQ(0, zero_bv->rank1(0));
    ASSERT_EQ(1, one_bv->rank1(0));
    size_t r = 0;
    for(size_t i = 0; i < random3_v.size(); i++)
    {
        r += random3_v[i];
        EXPECT_EQ(r, random3_bv->rank1(i)) << "when i = " << i;
    }
}

TEST_F(BitVectorTest, RankOutOfRange)
{
    ASSERT_THROW(zero_bv->rank1(5), std::out_of_range);
    ASSERT_THROW(random3_bv->rank1(random3_v.size()), std::out_of_range);
}

TEST_F(BitVectorTest, Select)
{
    ASSERT_EQ(0, zero_bv->select(0));
    ASSERT_EQ(1, one_bv->select(0));
    for(size_t i = 0; i < random3_v.size(); i++)
        EXPECT_EQ(random3_v[i], random3_bv->select(i)) << "when i = " << i;
}

TEST_F(BitVectorTest, SelectOutOfRange)
{
    ASSERT_THROW(zero_bv->select(5), std::out_of_range);
    ASSERT_THROW(random3_bv->select(random3_v.size()), std::out_of_range);
}

class WaveletTreeTest : public ::testing::Test
{
protected:
    WaveletTree *a_wt,
                *ab_wt,
                *abc_wt,
                *zero_wt,
                *test1_wt,
                *test2_wt,
                *test3_wt,
                *test4_wt,
                *test5_wt,
                *long_wt;
    std::string empty_str,
                a_str{"a"},
                ab_str{"ab"},
                abc_str{"abc"},
                zero_str{'\0'}, // Must use initialization list of chars when want to include '\0's.
                test1_str{"abcde"},
                test2_str{"abcdedcba"},
                test3_str{'\0', '1', '2', '3', '4', '5', 'a', 'b', 'c', 'd', 'd', 'd', 'd', 'd'},
                test4_str{'a', 'a', 'a', 'a', '\0'},
                test5_str{'a', 'a', 'b', 'b', '\0', 'd', 'd', 'e', 'e'},
                long_str{"One should not pursue goals that are easily achieved. "
                         "One must develop an instinct for \xFF what one can just barely "
                         "achieve through one's greatest efforts.\xE2"}; // NB: Extended ASCII

    virtual void SetUp(void)
    {
        a_wt = new WaveletTree(a_str, false);
        ab_wt = new WaveletTree(ab_str, false);
        abc_wt = new WaveletTree(abc_str, false);
        zero_wt = new WaveletTree(zero_str, false);
        test1_wt = new WaveletTree(test1_str, false);
        test2_wt = new WaveletTree(test2_str, false);
        test3_wt = new WaveletTree(test3_str, false);
        test4_wt = new WaveletTree(test4_str, false);
        test5_wt = new WaveletTree(test5_str, false);
        long_wt = new WaveletTree(long_str, false);
    }

    virtual void TearDown(void)
    {
        delete a_wt;
        delete ab_wt;
        delete abc_wt;
        delete zero_wt;
        delete test1_wt;
        delete test2_wt;
        delete test3_wt;
        delete test4_wt;
        delete test5_wt;
        delete long_wt;
    }

    static bool alphabet_matches(const std::string & alphabet,
                                 const std::string & s)
    {
        /* The idea of this function is to use logic that
         is likely to be somewhat independent of that
         which was used to generate the WaveletTree
         alphabet. Hence the odd, slow algorithm.
         */
        size_t total_matches = 0;
        for(std::string::const_iterator it = alphabet.begin(); it != alphabet.end(); it++)
        {
            size_t matches = 0;
            for(std::string::const_iterator jt = s.begin(); jt != s.end(); jt++)
                if(*jt == *it) matches++;
            if(matches == 0) return false;
            total_matches += matches;
        }
        return total_matches == s.size() ? true : false;
    }
};

TEST_F(WaveletTreeTest, Size)
{
    ASSERT_EQ(1, a_wt->size());
    ASSERT_EQ(2, ab_wt->size());
    ASSERT_EQ(3, abc_wt->size());
    ASSERT_EQ(1, zero_wt->size());
    ASSERT_EQ(test1_str.size(), test1_wt->size());
    ASSERT_EQ(test2_str.size(), test2_wt->size());
    ASSERT_EQ(test3_str.size(), test3_wt->size());
    ASSERT_EQ(test4_str.size(), test4_wt->size());
    ASSERT_EQ(test5_str.size(), test5_wt->size());
}

TEST_F(WaveletTreeTest, Alphabet)
{
    ASSERT_TRUE(alphabet_matches(a_wt->get_alphabet(), a_str));
    ASSERT_TRUE(alphabet_matches(ab_wt->get_alphabet(), ab_str));
    ASSERT_TRUE(alphabet_matches(abc_wt->get_alphabet(), abc_str));
    ASSERT_TRUE(alphabet_matches(zero_wt->get_alphabet(), zero_str));
    ASSERT_TRUE(alphabet_matches(test1_wt->get_alphabet(), test1_str));
    ASSERT_TRUE(alphabet_matches(test2_wt->get_alphabet(), test2_str));
    ASSERT_TRUE(alphabet_matches(test3_wt->get_alphabet(), test3_str));
    ASSERT_TRUE(alphabet_matches(test4_wt->get_alphabet(), test4_str));
    ASSERT_TRUE(alphabet_matches(test5_wt->get_alphabet(), test5_str));
    ASSERT_TRUE(alphabet_matches(long_wt->get_alphabet(), long_str));
}

TEST_F(WaveletTreeTest, CumFreq)
{
    std::string alphabet = long_wt->get_alphabet();
    for(std::string::iterator it = alphabet.begin(); it != alphabet.end(); it++)
    {
        size_t r = 0;
        for(std::string::iterator jt = long_str.begin(); jt != long_str.end(); jt++)
            if(static_cast<unsigned char>(*jt) < static_cast<unsigned char>(*it)) r++;
        EXPECT_EQ(r, long_wt->cum_freq(*it));
    }
}

TEST_F(WaveletTreeTest, Rank)
{
    ASSERT_EQ(1, a_wt->rank(0, 'a'));
    ASSERT_EQ(0, a_wt->rank(0, 'b'));
    ASSERT_EQ(1, zero_wt->rank(0, '\0'));
    ASSERT_EQ(0, zero_wt->rank(0, 'a'));
    for(int c = 0; c < 256; c++)
    {
        size_t r = 0;
        for(size_t i = 0; i < test5_wt->size(); i++)
        {
            if(test5_str[i] == (char) c) r++;
            EXPECT_EQ(r, test5_wt->rank(i, c)) << "when i = " << i << " and c = " << c;
        }
        r = 0;
        for(size_t i = 0; i < long_wt->size(); i++)
        {
            if(long_str[i] == (char) c) r++;
            EXPECT_EQ(r, long_wt->rank(i, c)) << "when i = " << i << " and c = " << c;
        }
    }
}

TEST_F(WaveletTreeTest, Empty)
{
    ASSERT_THROW(WaveletTree(empty_str, false), std::length_error);
}

TEST_F(WaveletTreeTest, RankOutOfRange)
{
    ASSERT_THROW(a_wt->rank(5, 'a'), std::out_of_range);
    ASSERT_THROW(long_wt->rank(long_wt->size(), 'a'), std::out_of_range);
}

TEST_F(WaveletTreeTest, Select)
{
    ASSERT_EQ('a', a_wt->select(0));
    ASSERT_EQ('\0', zero_wt->select(0));
    ASSERT_EQ('\0', test5_wt->select(4));
    for(size_t i = 0; i < long_str.size(); i++)
        EXPECT_EQ(long_str[i], long_wt->select(i)) << "when i = " << i;
}

TEST_F(WaveletTreeTest, SelectOutOfRange)
{
    ASSERT_THROW(a_wt->select(5), std::out_of_range);
    ASSERT_THROW(long_wt->select(long_wt->size()), std::out_of_range);
}

class FMIndexTest : public ::testing::Test
{
protected:
    FMIndex *zero_fmi,
            *aaaaa_fmi,
            *test_fmi,
            *long_fmi,
            *extra_fmi,
            *yet_another_fmi;
    const std::string zero_str{'\0'},
                      aaaaa_str{"aaaaa"},
                      test_str{'\0', 'a', 'b', 'c', 'd', 'e', '\0', 'h', 'e', 'l', 'l', 'o', '\xAB', // NB Extended ASCII
                               ' ', 't', 'h', 'e', 'r', 'e', ' ', 'a', 'n', 'd', ' ', 'h', 'e',
                               'l', 'l', 'o', ' ', 'g', 'o', 'o', 'd', 'b', 'y', 'e', 'h', 'e',
                               'l', 'l', 'o'}, // Must use initialization list of chars because of '\0's.
                      long_str{"Western civilization, it seems to me, stands by two great "
                               "heritages. One is the scientific spirit of adventure --- the "
                               "adventure into the unknown, an unknown which must be "
                               "recognized as being unknown in order to be explored; the "
                               "demand that the unanswerable mysteries of the universe remain "
                               "unanswered; the attitude that all is uncertain; to summarize it "
                               "--- the humility of the intellect. The other great heritage is "
                               "Christian ethics --- the basis of action on love, the brotherhood of "
                               "all men, the value of the individual --- the humility of the spirit."},
                      extra_str{"this\nshould\ncause\ntrouble"}, // Not really of course!
                      yet_another_str{"blah-de-blah"};
    std::list<std::pair<FMIndex::const_iterator, FMIndex::const_reverse_iterator>> matches;

    virtual void SetUp(void)
    {
        zero_fmi = new FMIndex(zero_str);
        aaaaa_fmi = new FMIndex(aaaaa_str);
        test_fmi = new FMIndex(test_str);
        long_fmi = new FMIndex(long_str);
        extra_fmi = new FMIndex(extra_str);
        yet_another_fmi = new FMIndex(yet_another_str);
    }

    virtual void TearDown(void)
    {
        delete zero_fmi;
        delete aaaaa_fmi;
        delete test_fmi;
        delete long_fmi;
        delete extra_fmi;
        delete yet_another_fmi;
    }

    static std::string scan_back(FMIndex::const_reverse_iterator & it,
                                 const size_t len)
    {
        std::ostringstream ss;
        std::copy_n(it, len, std::ostream_iterator<char>(ss));
        std::string s = ss.str();
        return std::string(s.rbegin(), s.rend());
    }

    static std::string scan_forward(FMIndex::const_iterator & it,
                                    const size_t len)
    {
        std::ostringstream ss;
        std::copy_n(it, len, std::ostream_iterator<char>(ss));
        return ss.str();
    }

    static std::string get_text(const FMIndex & fmi)
    {
        std::ostringstream ss;
        std::copy_n(fmi.begin(), fmi.size(), std::ostream_iterator<char>(ss));
        return ss.str();
    }
};

TEST_F(FMIndexTest, Empty)
{
    ASSERT_THROW(FMIndex{std::string()}, std::length_error);
    ASSERT_THROW(long_fmi->find(matches, std::string()), std::length_error);
}

TEST_F(FMIndexTest, Find)
{
    ASSERT_EQ(1, zero_fmi->find(matches, std::string{'\0'}));
    ASSERT_EQ(5, aaaaa_fmi->find(matches, std::string("a")));
    ASSERT_EQ(3, test_fmi->find(matches, std::string("hello")));
    ASSERT_EQ(17, long_fmi->find(matches, std::string("the")));
    ASSERT_EQ(1, long_fmi->find_lines(std::string("Western civilization, it")).size());
    ASSERT_EQ(2, extra_fmi->find_lines(std::string("t")).size());
    ASSERT_EQ(1, yet_another_fmi->findn(std::string("-de")));
}

TEST_F(FMIndexTest, Scan)
{
    long_fmi->find(matches, std::string("individual"));
    ASSERT_EQ(1, matches.size());
    ASSERT_EQ("value of the ", scan_back(matches.front().second, 13));
    ASSERT_EQ(" --- the humility", scan_forward(matches.front().first, 17));
    ASSERT_EQ(long_str, get_text(*long_fmi));
}

TEST(Serializing, Basic)
{
    const uint64_t x = 13446544033719551615LLU; // Random value.
    std::ostringstream s;
    serialize_as_chars(std::ostreambuf_iterator<char>(s), x);
    uint64_t y;
    std::istringstream ss(s.str());
    deserialize_from_chars(std::istreambuf_iterator<char>(ss), y);
    ASSERT_EQ(x, y);
}

TEST_F(BitVectorTest, Serializing)
{
    std::ostringstream s;
    random3_bv->serialize(std::ostreambuf_iterator<char>(s));
    std::istringstream ss(s.str());
    BitVector bv{std::istreambuf_iterator<char>(ss)}; // Avoid "most vexing parse"
    ASSERT_EQ(random3_bv->size(), bv.size());
    for(size_t i = 0; i < bv.size(); i++)
        EXPECT_EQ(random3_bv->select(i), bv.select(i));
    for(size_t i = 0; i < bv.size(); i++)
        EXPECT_EQ(random3_bv->rank1(i), bv.rank1(i));
}

TEST_F(WaveletTreeTest, Serializing)
{
    std::ostringstream s;
    long_wt->serialize(std::ostreambuf_iterator<char>(s));
    std::istringstream ss(s.str());
    WaveletTree wt{std::istreambuf_iterator<char>(ss)}; // Avoid "most vexing parse"
    ASSERT_EQ(long_wt->size(), wt.size());
    for(size_t i = 0; i < wt.size(); i++)
        EXPECT_EQ(long_wt->select(i), wt.select(i));
}

TEST_F(FMIndexTest, Serializing)
{
    std::ostringstream s;
    long_fmi->serialize(std::ostreambuf_iterator<char>(s));
    std::istringstream ss(s.str());
    FMIndex fmi{std::istreambuf_iterator<char>(ss)}; // Avoid "most vexing parse"
    ASSERT_EQ(long_fmi->size(), fmi.size());
    ASSERT_EQ(long_str, get_text(fmi));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
