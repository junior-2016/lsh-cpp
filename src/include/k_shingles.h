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

        K_shingling(const ValueType &_value, const WeightType &_weight) : _value(_value), _weight(_weight) {}

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
            (*(result.emplace(substr, 0)).first)._weight++;
        }
        return result;
    }

    /**
     * 为处理DNA数据压缩的shingling结构体.
     * A T C G 用 two-bit 编码,内部结构体用 std::bitset
     * @tparam k k-shingling长度
     * @tparam flag 标识是否拥有权重数据,默认无权重(no_weight)
     */
    enum class WeightFlag : uint32_t {
        has_weight = 1,
        no_weight = 2
    };

    /**
     * DNA_Shingling 模板
     * @tparam k k_shingling length.
     * @tparam flag 只允许 WeighFlag::has_weight / WeightFlag::no_weight 两个值,不允许其他类型的值.
     * 注意 :
     * 虽然 WeightFlag 的底层类型是 uint32_t, 但是你用 uint32_t 来实例化也是不行的, 比如 DNA_Shingling<7,(uint32_t)(1)> 不能通过编译.
     * 这样就确保了强类型安全.
     */
    template<size_t k, auto flag>
    struct DNA_Shingling {
        static_assert(std::is_same_v<decltype(flag), WeightFlag>);
        static_assert((flag == WeightFlag::no_weight) || (flag == WeightFlag::has_weight));
    };

    template<size_t k>
    struct DNA_Shingling<k, WeightFlag::no_weight> {
        static constexpr size_t bits = 2 * k; // dna k-shingling 每一个位置用2-bit编码
        using ValueType = std::bitset<bits>;

        ValueType _value;

        explicit DNA_Shingling(const ValueType &_value) : _value(_value) {}

        [[nodiscard]] ValueType value() const { return _value; }

        bool operator==(const DNA_Shingling<k, WeightFlag::no_weight> &shingling) const {
            return _value == shingling._value;
        }
    };

    template<size_t k>
    struct DNA_Shingling<k, WeightFlag::has_weight> {
        static constexpr size_t bits = 2 * k; // dna k-shingling 每一个位置用2-bit编码
        using WeightType = uint32_t;
        using ValueType = std::bitset<bits>;

        mutable WeightType _weight;
        ValueType _value;

        DNA_Shingling(const ValueType &_value, const WeightType &_weight) : _value(_value), _weight(_weight) {}

        [[nodiscard]] ValueType value() const { return _value; }

        [[nodiscard]] WeightType weight() const { return _weight; }

        bool operator==(const DNA_Shingling<k, WeightFlag::has_weight> &shingling) const {
            return _value == shingling._value;
        }
    };

    // dna_shingling 编码压缩
    template<size_t k, auto flag>
    typename DNA_Shingling<k, flag>::ValueType dna_shingling_encode(const std::string_view &dna_shingling) {
        using ValueType = typename DNA_Shingling<k, flag>::ValueType;

        // std::bitset 初始时全0, bitset::set(pos,true/false)中pos是从右边最低位开始,
        // 比如 bitset<5> temp; temp.set(0,true); => temp.to_string() == "00001"
        ValueType value;

        size_t pos = DNA_Shingling<k, flag>::bits - 1; // 从最高位开始
        for (const auto &ch: dna_shingling) {
            switch (ch) {
                case 'A': // 00,不用设置1,但需要修改pos
                    pos -= 2;
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
                default:
                    break;
            }
        }
        return value;
    }

    // 注意,如果源字符串长度小于k,比如 str = AAAT, k = 5, 那么编码的时候依然会按 k = 5 编码为 00_00_00_01_00, 解码的时候就变成 AAATA.
    // 要解决这个问题需要多引入一个member表示字符串的真正长度,但考虑到极大多数的dna都非常长,并且长度远大于k,所以这里不考虑这个问题.
    template<size_t k, auto flag>
    std::string dna_shingling_decode(const DNA_Shingling<k, flag> &shingling) {
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

    template<size_t k, auto flag>
    HashSet <DNA_Shingling<k, flag>> split_dna_shingling(const std::string_view &string) {
        if (k >= string.size()) { // 这里的string.size()虽然是constexpr,但只有当string是编译期确定才有效.所以不能用if constexpr.
            if constexpr (flag == WeightFlag::has_weight) {
                return {DNA_Shingling<k, flag>{dna_shingling_encode<k, flag>(string), 1}};
            } else {
                return {DNA_Shingling<k, flag>{dna_shingling_encode<k, flag>(string)}};
            }
        }
        size_t N = string.size() - k + 1;
        HashSet <DNA_Shingling<k, flag>> result(N); // 预先设置哈希表容器大小避免扩容开销
        for (size_t i = 0; i < N; i++) {
            auto substr = string.substr(i, k);
            if constexpr (flag == WeightFlag::has_weight) {
                // 使用 emplace 在容器里就地构造对象,避免构造后拷贝开销.
                (*(result.emplace(dna_shingling_encode<k, flag>(substr), 0)).first)._weight++;
            } else {
                result.emplace(dna_shingling_encode<k, flag>(substr));
            }
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

    // dna shingling 用std::hash(bitset)得到哈希值.
    // TODO: 理论上也可以直接 bitset cast to size_t 作为哈希值,这样更快.
    //  另外 bitset cast 为 size_t 还可以直接作为权重最小哈希向量(Weighted MinHash Vector)的位置编码..
    template<size_t k, auto flag>
    struct hash<LSH_CPP::DNA_Shingling<k, flag>> {
        std::size_t operator()(LSH_CPP::DNA_Shingling<k, flag> const &dna_shingling) const {
            return phmap::HashState::combine(0, dna_shingling._value);
        }
    };
}
#endif //LSH_CPP_K_SHINGLES_H
