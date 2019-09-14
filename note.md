## C++ programming 
### std::vector initialize
- 如果提前知道vector需要多少内存且无需初始化,并只用做一次值更新,
最优的策略是使用 reserve(N) 直接分配内存且不初始化,
后面连续用 push_back(v) 分配你需要的值,
push_back不会做任何扩容操作,节省了时间.
需要注意的是,如果用了reserve(N),分配值的时候必须用push_back,
不能用operator[]. 
- 如果vector需要初始化,初始化的值将用在后面更新值的比较中; 
或者vector更新值不是顺序的; 
或者vector可能多次调用函数分配赋值而不是仅赋值一次, 
这三种情况就不能采用 reserve + push_back的方案, 
而应该用 resize + operator[] 的方法.
- 将一个计算好的vector传递给结构体内的vector成员,
通过vector_member.swap(vector_compute)是最快的,只需要常数时间,
因为这个操作仅仅交换两个vector的data_pointer而已.
```c++
struct Test{
   vector<int> member;
   void update(){
       vector<int> compute;
       // compute ...
       member.swap(compute);
   }
}
```

### C++17 parallelism TS
gcc 9.1.0 基于 intel TBB 实现了C++17的并行策略提案, 对于常见的std::sort,
std::for_each,std::reduce...等算法,可以提供 execution_tag(seq/par/par_unseq)
来实现并行计算. 需要注意的是几种 execution_tag 的使用场景,要尽量避免data race.
<br>
reference:<br>
(1) [https://devblogs.microsoft.com/cppblog/using-c17-parallel-algorithms-for-better-performance/]<br>
(2) [https://www.bfilipek.com/2017/08/cpp17-details-parallel.html]<br>
(3) [https://www.bfilipek.com/2018/11/parallel-alg-perf.html]<br>
(4) [https://en.cppreference.com/w/cpp/experimental/parallelism]<br>
(5) [https://en.cppreference.com/w/cpp/algorithm/execution_policy_tag_t]<br>
(6) [https://stackoverflow.com/questions/39954678/difference-between-execution-policies-and-when-to-use-them]<br>

### git notes
push local branch to remote branch with different name:
```bash
$ git push -u origin my_branch:remote_branch
```

### STL container insert & emplace
- 如果是就地插入(右值构造一个对象并且立即插入),应该考虑用 emplace 直接在容器内找一个
正确的位置申请空间并创建对象,从而避免预先创建对象然后拷贝/移动到容器中这些繁琐的动作.
如果插入的元素是容器中已经存在的,不会重复构造一次.注意使用 set<A>.emplace(args..) 传递参数时,结构体A必须有相对应的构造函数,并且
传递的参数顺序必须一致,不要传递错误. 
- 对于 tree_set/tree_map(c++普通的set/map),可以通过 itr = lower_bound(search_key);
 if(itr->first!=search_key) m.emplace_hint(itr,args...) 
 来避免反复搜索,因为前面 lower_bound 已经搜索过一次了,利用好这个信息可以加速插入.
 但是哈希set/哈希map不能这样做,因为哈希方法不是平衡树那样的顺序容器,
 是通过hash value来确定位置的,如果查不到就是map.end()了,没有有效的信息可以利用,
 也没有lower_bound/upper_bound接口.
- 对于 map<key,value> 的 emplace 构造,有两种方式,一种是:
 map<key,value>.emplace(key_construct{args},value_construct{args})
 这种方法的好处是简洁易懂,但是坏处也很明显:只是一部分实现了容器内 in-place construct,
 因为必须我们预先构造key和value(这个构造并不是容器内构造),然后构造 std::pair<key,value> 
 (这个对象的构造才是容器里就地构造的).如何实现key,value对象本身也是在容器里构造呢,这就需要借助
 std::piecewise_construct_t(c++17里面可以用try_emplace来代替,但是语法上不够简洁易懂)
 所以第二种方法是: map.emplace(std:: piecewise_construct,
 std:: forward_as_tuple(key_args),std:: forward_as_tuple(value_args));
 或者 map.try_emplace(key_args,... , value_args...).
- reference:
<br>[http://clarkkromenaker.com/post/cpp-emplace/]
<br>[https://stackoverflow.com/questions/26446352/what-is-the-difference-between-unordered-map-emplace-and-unordered-map-ins]
<br>[https://en.cppreference.com/w/cpp/container/unordered_map/emplace]
<br>[https://en.cppreference.com/w/cpp/container/map/try_emplace]
<br>[http://www.cplusplus.com/reference/unordered_set/unordered_set/emplace/]
<br>[https://stackoverflow.com/questions/41507671/what-is-the-use-of-emplace-hint-in-map]
