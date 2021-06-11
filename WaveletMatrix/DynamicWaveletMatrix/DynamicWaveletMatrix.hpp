#include <bits/stdc++.h>
#include "../DynamicBitVector/DynamicBitVector.hpp"

class DynamicWaveletMatrix {
public:
    std::vector<DynamicBitVector> bit_arrays;
    std::vector<uint64_t> begin_one;               // 各bitの1の開始位置

    uint64_t size;                                 // 与えられた配列のサイズ
    uint64_t maximum_element;                      // 最大の文字
    uint64_t bit_size;                             // 文字を表すのに必要なbit数

public:
    // max_element: 入ってくる中で一番大きい数値
    DynamicWaveletMatrix (uint64_t maximum_element) : size(0), maximum_element(maximum_element + 1) {
        this->bit_size = this->get_num_of_bit(maximum_element);
        if (bit_size == 0) {
            bit_size = 1;
        }
        this->begin_one.resize(bit_size);

        for (uint64_t i = 0; i < bit_size; ++i) {
            DynamicBitVector sv;
            bit_arrays.push_back(sv);
        }
    }

    DynamicWaveletMatrix (uint64_t num_of_alphabet, const std::vector<uint64_t> &array) : size(0), maximum_element(num_of_alphabet + 1) {
        this->bit_size = this->get_num_of_bit(num_of_alphabet);
        if (bit_size == 0) {
            bit_size = 1;
        }
        this->begin_one.resize(bit_size);

        if (array.empty()) {
            for (uint64_t i = 0; i < bit_size; ++i) {
                DynamicBitVector sv;
                bit_arrays.push_back(sv);
            }
            return;
        }

        size = array.size();

        std::vector<uint64_t> v(array), b(array.size(), 0);

        for (uint64_t i = 0; i < bit_size; ++i) {

            std::vector<uint64_t> temp;
            // 0をtempにいれてく
            for (uint64_t j = 0; j < v.size(); ++j) {
                uint64_t c = v.at(j);
                uint64_t bit = (c >> (bit_size - i - 1)) & 1;  //　上からi番目のbit
                if (bit == 0) {
                    temp.push_back(c);
                    b[j] = 0;
                }
            }

            this->begin_one.at(i) = temp.size();

            // 1をtempにいれてく
            for (uint64_t j = 0; j < v.size(); ++j) {
                uint64_t c = v.at(j);
                uint64_t bit = (c >> (bit_size - i - 1)) & 1;  //　上からi番目のbit
                if (bit == 1) {
                    temp.push_back(c);
                    b[j] = 1;
                }
            }

            DynamicBitVector dbv(b);
            bit_arrays.emplace_back(dbv);
            v = temp;
        }
    }

    // v[pos]
    uint64_t access(uint64_t pos) {
        if (pos >= this->size) { return NOTFOUND; }

        uint64_t c = 0;
        for (uint64_t i = 0; i < bit_arrays.size(); ++i) {
            uint64_t bit = bit_arrays.at(i).access(pos);   // もとの数値がのi番目のbit
            c = (c <<= 1) | bit;
            pos = bit_arrays.at(i).rank(bit, pos);
            if (bit) {
                pos += this->begin_one.at(i);
            }
        }
        return c;
    }

    // v[0, pos)のcの数
    uint64_t rank(uint64_t c, uint64_t pos) {
        assert(pos <= size);
        if (c >= maximum_element) {
            return 0;
        }

        uint64_t left = 0, right = pos;
        for (uint64_t i = 0; i < bit_size; ++i) {
            const uint64_t bit = (c >> (bit_size - i - 1)) & 1;  // 上からi番目のbit
            left = bit_arrays.at(i).rank(bit, left);             // cのi番目のbitと同じ数値の数
            right = bit_arrays.at(i).rank(bit, right);           // cのi番目のbitと同じ数値の数
            if (bit) {
                left += this->begin_one.at(i);
                right += this->begin_one.at(i);
            }
        }

        return right - left;
    }

    // i番目のcの位置 + 1を返す。rankは1-origin
    uint64_t select(uint64_t c, uint64_t rank) {
        assert(rank > 0);
        if (c >= maximum_element) {
            return NOTFOUND;
        }

        uint64_t left = 0;
        for (uint64_t i = 0; i < bit_size; ++i) {
            const uint64_t bit = (c >> (bit_size - i - 1)) & 1;  // 上からi番目のbit
            left = bit_arrays.at(i).rank(bit, left);               // cのi番目のbitと同じ数値の数
            if (bit) {
                left += this->begin_one.at(i);
            }
        }

        uint64_t index = left + rank;
        for (uint64_t i = 0; i < bit_arrays.size(); ++i){
            uint64_t bit = ((c >> i) & 1);      // 下からi番目のbit
            if (bit == 1) {
                index -= this->begin_one.at(bit_size - i - 1);
            }
            index = this->bit_arrays.at(bit_size - i - 1).select(bit, index);
        }
        return index;
    }

    // posにcを挿入する
    void insert(uint64_t pos, uint64_t c) {
        DCHECK_LE(pos, this->size);

        for (uint64_t i = 0; i < bit_arrays.size(); ++i) {
            const uint64_t bit = (c >> (bit_size - i - 1)) & 1;  //　上からi番目のbit
            bit_arrays.at(i).insert(pos, bit);
            pos = bit_arrays.at(i).rank(bit, pos);
            if (bit) {
                pos += this->begin_one.at(i);
            }
            else {
                this->begin_one.at(i)++;
            }
        }

        this->size++;
    }

    // 末尾にcを追加する
    void push_back(uint64_t c) {
        this->insert(this->size, c);
    }

    // posを削除する
    void erase(uint64_t pos) {
        assert(pos < this->size);
        if (pos >= this->size) {
            throw "Segmentation fault";
        }

        for (uint64_t i = 0; i < bit_arrays.size(); ++i) {
            uint64_t bit = bit_arrays.at(i).access(pos);   // もとの数値のi番目のbit

            auto next_pos = bit_arrays.at(i).rank(bit, pos);
            bit_arrays.at(i).erase(pos);

            if (bit) {
                next_pos += this->begin_one.at(i);
            }
            else {
                this->begin_one.at(i)--;
            }
            pos = next_pos;
        }
        this->size--;
    }

    // posにcをセットする
    void update(uint64_t pos, uint64_t c) {
        assert(pos < this->size);
        this->erase(pos);
        this->insert(pos, c);
    }

    // 他の操作は通常のWavelet Matrixと同じ

    // T[begin_pos, end_pos)でx <= c < yを満たす最大のcを返す
    uint64_t prevValue(const uint64_t begin_pos, const uint64_t end_pos, const uint64_t x, uint64_t y) {
        assert(end_pos <= size);
        const uint64_t num = end_pos - begin_pos;

        if (x >= y or y == 0) {
            return NOTFOUND;
        }
        if (y > maximum_element) {
            y = maximum_element;
        }

        if (begin_pos >= end_pos) {
            return NOTFOUND;
        }
        if (x >= maximum_element || end_pos == 0) {
            return NOTFOUND;
        }

        y--; // x <= c <= yにする

        std::stack<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, bool>> s;   // (begin, end, depth, c, tight)
        s.emplace(std::make_tuple(begin_pos, end_pos, 0, 0, true));

        while (not s.empty()) {
            uint64_t b, e, depth, c;
            bool tight;
            std::tie(b, e, depth, c, tight) = s.top(); s.pop();

            if (depth == bit_size) {
                if (c >= x) {
                    return c;
                }
                continue;
            }

            const uint64_t bit = (y >> (bit_size - depth - 1)) & 1;

            const uint64_t rank0_begin = this->bit_arrays.at(depth).rank(0, b);
            const uint64_t rank0_end = this->bit_arrays.at(depth).rank(0, e);
            const uint64_t rank1_begin = b - rank0_begin;
            const uint64_t rank1_end = e - rank0_end;

            // d番目のbitが0のものを使う
            const uint64_t b0 = rank0_begin;
            const uint64_t e0 = rank0_end;
            if (b0 != e0) { // 範囲がつぶれてない
                const uint64_t c0 = ((c << 1) | 0);
                s.emplace(std::make_tuple(b0, e0, depth + 1, c0, tight and bit == 0));
            }

            // d番目のbitが1のものを使う
            const uint64_t b1 = this->begin_one.at(depth) + rank1_begin;
            const uint64_t e1 = this->begin_one.at(depth) + rank1_end;
            if (b1 != e1) {
                if (not tight or bit == 1) {
                    const auto c1 = ((c << 1) | 1);
                    s.emplace(std::make_tuple(b1, e1, depth + 1, c1, tight));
                }
            }
        }

        return NOTFOUND;
    }

    // T[begin_pos, end_pos)でx <= c < yを満たす最小のcを返す
    uint64_t nextValue(const uint64_t begin_pos, const uint64_t end_pos, const uint64_t x, const uint64_t y) {
        assert(end_pos <= size);
        const uint64_t num = end_pos - begin_pos;

        if (x >= y or y == 0) {
            return NOTFOUND;
        }

        if (begin_pos >= end_pos) {
            return NOTFOUND;
        }
        if (x >= maximum_element || end_pos == 0) {
            return NOTFOUND;
        }

        std::stack<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t, bool>> s;   // (begin, end, depth, c, tight)
        s.emplace(std::make_tuple(begin_pos, end_pos, 0, 0, true));

        while (not s.empty()) {
            uint64_t b, e, depth, c;
            bool tight;
            std::tie(b, e, depth, c, tight) = s.top(); s.pop();

            if (depth == bit_size) {
                if (c < y) {
                    return c;
                }
                continue;
            }

            const uint64_t bit = (x >> (bit_size - depth - 1)) & 1;

            const uint64_t rank0_begin = this->bit_arrays.at(depth).rank(0, b);
            const uint64_t rank0_end = this->bit_arrays.at(depth).rank(0, e);
            const uint64_t rank1_begin = b - rank0_begin;
            const uint64_t rank1_end = e - rank0_end;

            // d番目のbitが1のものを使う
            const uint64_t b1 = this->begin_one.at(depth) + rank1_begin;
            const uint64_t e1 = this->begin_one.at(depth) + rank1_end;
            if (b1 != e1) {
                const auto c1 = ((c << 1) | 1);
                s.emplace(std::make_tuple(b1, e1, depth + 1, c1, tight and bit == 1));
            }

            // d番目のbitが0のものを使う
            const uint64_t b0 = rank0_begin;
            const uint64_t e0 = rank0_end;
            if (b0 != e0) {
                if (not tight or bit == 0) {
                    const uint64_t c0 = ((c << 1) | 0);
                    s.emplace(std::make_tuple(b0, e0, depth + 1, c0, tight));
                }
            }
        }

        return NOTFOUND;
    }


private:
    uint64_t get_num_of_bit(uint64_t x) {
        if (x == 0) return 0;
        x--;
        uint64_t bit_num = 0;
        while (x >> bit_num) {
            ++bit_num;
        }
        return bit_num;
    }
};
