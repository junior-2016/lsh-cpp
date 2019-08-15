//
// Created by junior on 2019/8/10.
//

#ifndef LSH_CPP_K_SHINGLES_H
#define LSH_CPP_K_SHINGLES_H

#include "lsh_cpp.h"
#include "util.h"

namespace LSH_CPP {
    struct K_shingling {
        using WeightType = uint32_t; //size_t;
        using ValueType = std::string_view;

        ValueType _value;
        mutable WeightType _weight;

        // TODO: pos_list 后面可以考虑用BitSet压缩,不然内存开销会很大. pos_list
        //  方法: 先使用相对位置编码(一阶差分),保留pos_list[0],剩下的保存相邻元素的差值,比如原pos_list是{ 0,10,28,35,45 }
        //  相对位置编码为 { 0,10-0,28-10,35-28,45-35 } = { 0,10,18,7,10 } .
        //  然后再用unary_bit_set对上面的序列做位编码(见<<信息检索>>相关压缩方法).
        //  2019/8/12 pos_list 已经严重影响内存开销,先考虑移除(后面再做压缩)...
        // mutable std::vector<size_t> pos_list;

        // 通过weight()方法得到权重(重复元素个数).
        WeightType weight() const {
            // return pos_list.size();
            return _weight;
        }

        // 通过value()方法获取内部元素(k_shingling 串)
        ValueType value() const {
            return _value;
        }

        bool operator==(const K_shingling &k_shingling) const {
            return _value == k_shingling._value;
        }
    };

    /**
     *
     * split string to k_mer.
     * @param string 源字符串,考虑用 string_view 优化,因为string_view只储存一个data pointer和一个字符串长度,
     * 这样使用substr()将不会有任何的内存开销,速度也相对更快,只需要保证源字符串自始至终都存在(不被析构)就好.
     * 参考: https://stackoverflow.com/questions/42921998/c-efficiently-get-substring-of-string-with-index
     * 经过比较,使用std::string_view比直接用std::string快了一倍!
     * 但是string_view毕竟只是持有原字符串的地址,所以一旦发生源字符串的析构就会导致严重错误,比如我下面的调用代码:
     * int main(){
     *     std::vector<std::string_view> ret;
     *     {
     *        std::string s = "abcdefghijklmnopqrstuvwxyz";
     *        ret = split_k_mer<32>(s, 3);
     *     } // 离开作用域后源字符串s被析构
     *     for (const auto& item:ret){
     *        std::cout<<item<<" ";
     *        // 输出ret的元素,这时string_view持有的地址所在的内容已经析构或者被别的方式处理了,
     *        // 此时输出结果是 undefined behavior (未定义行为).
     *     }
     * }
     * 因此在使用string_view时要格外小心源字符串的生命周期
     */
    // 更快的k_mer_split,但是需要考虑源字符串的生命周期.
    HashSet <K_shingling> split_k_shingling_fast(const std::string_view &string, size_t k) {
        if (k >= string.size()) {
            // return {{string, {0}}};
            return {{string, 1}};
        }
        size_t N = string.size() - k + 1;
        HashSet<K_shingling> result;
        for (size_t i = 0; i < N; i++) {
            auto substr = string.substr(i, k);
            // 这里pos_list一开始留空,等后面插入(不管成功还是失败)后,从返回值的first得到迭代器,然后再修改pos_list.
            // 本来 unordered_set.insert(Key) 返回的应该是 { const_iterator, bool },
            // 但因为我们把 K_mer 结构体的 pos_list 成员改为 mutable(可变). 因此,即使是const_iterator,也可以修改它.
            // 这样写可以最大效率实现插入和修改pos_list,完全不需要插入一次,再查找一次,因为插入的时候本身就是在查找.
            // 另外,这种写法完全可以用 map [key(不可变部分)] = value (可变部分) 来代替.
            // (*(result.insert(K_shingling{substr, {}})).first).pos_list.push_back(i);
            (*(result.insert(K_shingling{substr, 0})).first)._weight++;
        }
        return result;
    }

    /**
     * 为处理DNA数据压缩的shingling结构体.
     * A T C G 用 two-bit 编码,内部结构体用 std::bitset
     * @tparam k k-shingling长度
     */
    template<size_t k>
    struct DNA_Shingling {
        static constexpr size_t bits = 2 * k; // dna k-shingling 每一个位置用2-bit编码
        // using WeightType = uint32_t;
        using ValueType = std::bitset<bits>;

        // mutable WeightType _weight;
        ValueType _value;

//        DNA_Shingling(const ValueType &_value, const WeightType &_weight) : _value(_value), _weight(_weight) {}
//
//        [[nodiscard]] WeightType weight() const { return _weight; }

        [[nodiscard]] ValueType value() const { return _value; }

        bool operator==(const DNA_Shingling<k> &shingling) const {
            return _value == shingling._value;
        }
    };

    // dna_shingling 编码压缩
    template<size_t k>
    typename DNA_Shingling<k>::ValueType dna_shingling_encode(const std::string_view &dna_shingling) {
        using ValueType = typename DNA_Shingling<k>::ValueType;

        // std::bitset 初始时全0, bitset::set(pos,true/false)中pos是从右边最低位开始,
        // 比如 bitset<5> temp; temp.set(0,true); => temp.to_string() == "00001"
        ValueType value;

        size_t pos = DNA_Shingling<k>::bits - 1; // 从最高位开始
        for (const auto &ch: dna_shingling) {
            switch (ch) {
                case 'A': // 00,不用设置1,但需要修改pos
                    pos--;
                    pos--;
                    break;
                case 'T': // 01
                    pos--;
                    value.set(pos--, 1);
                    break;
                case 'C': // 10
                    value.set(pos--, 1);
                    pos--;
                    break;
                case 'G': // 11
                    value.set(pos--, 1);
                    value.set(pos--, 1);
                    break;
                default: // TODO 如果不在 A/T/C/G 四个字符范围,可能需要考虑抛出异常(但是会影响性能)
                    break;
            }
        }
        return value;
    }

    // 注意,如果源字符串长度小于k,比如 str = AAAT, k = 5, 那么编码的时候依然会按 k = 5 编码为 00_00_00_01_00, 解码的时候就变成 AAATA.
    // 要解决这个问题需要多引入一个member表示字符串的真正长度,但考虑到极大多数的dna都非常长,并且长度远大于k,所以这里不考虑这个问题.
    template<size_t k>
    std::string dna_shingling_decode(const DNA_Shingling<k> &shingling) {
        std::string dna;
        int pos = shingling.bits - 1; // 注意这里pos用int类型,否则下面的while会死循环.
        auto value = shingling.value();
        while (pos >= 0) {
            if (value[pos] == 0) {
                if (value[pos - 1] == 0) dna += 'A'; else dna += 'T';
            } else {
                if (value[pos - 1] == 0) dna += 'C'; else dna += 'G';
            }
            pos -= 2;
        }
        return dna;
    }

    template<size_t k>
    HashSet <DNA_Shingling<k>> split_dna_shingling(const std::string_view &string) {
        if (k >= string.size()) {
//            return {DNA_Shingling<k>{dna_shingling_encode<k>(string), 1}};
            return {DNA_Shingling<k>{dna_shingling_encode<k>(string)}};
        }
        size_t N = string.size() - k + 1;
        HashSet<DNA_Shingling<k>> result;
        for (size_t i = 0; i < N; i++) {
            auto substr = string.substr(i, k);
//            (*(result.insert(DNA_Shingling<k>{dna_shingling_encode<k>(substr), 0})).first)._weight++;
            result.insert(DNA_Shingling<k>{dna_shingling_encode<k>(substr)});
        }
        return result;
    }
}
namespace std {
    // inject specialization of std::hash for K_shingling
    template<>
    struct hash<LSH_CPP::K_shingling> {
        std::size_t operator()(LSH_CPP::K_shingling const &k_shingling) const {
            return phmap::HashState::combine(0, k_shingling._value);
        }
    };

    template<size_t k>
    struct hash<LSH_CPP::DNA_Shingling<k>> {
        std::size_t operator()(LSH_CPP::DNA_Shingling<k> const &dna_shingling) const {
            return phmap::HashState::combine(0, dna_shingling._value);
        }
    };
}
#endif //LSH_CPP_K_SHINGLES_H
