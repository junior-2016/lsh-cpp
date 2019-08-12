//
// Created by junior on 2019/8/12.
//

#ifndef LSH_CPP_LRU_CACHE_H
#define LSH_CPP_LRU_CACHE_H

#include "lsh_cpp.h"
#include "util.h"

namespace LSH_CPP {
    // reference: https://github.com/lamerman/cpp-lru-cache
    template<
            typename K,
            typename V,
            typename Hash,
            typename Eq,
            typename Alloc,
            template<typename/*K*/, typename/*V*/, typename /*Hash*/, typename /*Eq*/, typename /*Alloc*/> typename Map
    >
    class lru_cache {
    private:
        using k_v_pair_t = std::pair<K, V>;
        using Cache_List_t = std::list<k_v_pair_t>;
        using Cache_List_Iterator_t = typename Cache_List_t::iterator;
        using Cache_Map_t = Map<K, Cache_List_Iterator_t, Hash, Eq, Alloc>;
        using Cache_Map_Iterator_t = typename Cache_Map_t::iterator;
        size_t _max_size;
        Cache_Map_t cache_map;
        Cache_List_t cache_list;
    public:
        explicit lru_cache(size_t _max_size) : _max_size(_max_size) {}

        [[nodiscard]] size_t max_size() const { return _max_size; }

        bool contains(const K &key) const { return cache_map.find(key) != cache_map.end(); }

        void put(const K &key, const V &value) {
            cache_list.push_front({key, value});
            if (Cache_Map_Iterator_t it = cache_map.find(key);
                    it != cache_map.end()) {
                cache_list.erase(it->second);
                cache_map.erase(it); // erase by map::iterator
            }
            cache_map[key] = cache_list.begin();
            if (cache_map.size() > _max_size) { // 超过最大元素个数,移除list最后一个元素,同步删除map里的元素
                Cache_List_Iterator_t last = cache_list.end();
                last--; // 得到最后一个元素的迭代器
                cache_map.erase(last->first); // erase by key
                cache_list.pop_back();
            }
        }

        /**
         * use case:
         * if (auto ret = cache.get(key); ret.has_value()){
         *     // process *ret
         * } else {
         *    // ret is null option. other process ...
         * }
         * ps: 上面的use-case可以避免多次检查key.
         * 如果你用 if (cache.contains(key)) { auto ret = cache.get(key); }
         * contains 和 get 方法会对 key 检查两次,多了一点开销
         */
        std::optional<V> get(const K &key) {
            if (Cache_Map_Iterator_t it = cache_map.find(key);
                    it == cache_map.end()) {
                return std::nullopt; // 返回空值
            } else {
                // 将get()得到的元素移到list的队首.
                cache_list.splice(cache_list.begin(), cache_list, it->second);
                return it->second->second;
            }
        }
    };
}
#endif //LSH_CPP_LRU_CACHE_H
