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