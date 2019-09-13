//
// Created by junior on 2019/9/8.
//

#ifndef LSH_CPP_LSH_COSINE_SIMILARITY_H
#define LSH_CPP_LSH_COSINE_SIMILARITY_H
// 基于余弦相似度的 LSH 加速.
// TODO: K-mer natural vector and its application to the phylogenetic analysis of genetic sequences
//      构造 k-mer natural vector, 然后基于 cos similarity 计算相似度, 用 lsh 加速.
//      用这种方法与 weight minhash 比较结果.
// TODO
//      weight minhash 与 k-mer natural vector 都需要一个初始的全域向量(记录k_mer在sequence里出现的次数),
//      而这个过程需要得到每一个k-mer的位置编码.
//      之前实现的通用的 weight-minhash , 一种是在小dimension里直接提供权重向量,位置编码当然也是用户决定;
//      一种是大dimension里内部通过遍历次序确定权重向量和位置.
//      这里我后续想修改的是: 不再通过 dimension大小 区分实现,而是提供两种实现,不同实现用不同名称,
//      因为若dimension实在过大就理应使用第二个实现,不需要模板内部决定,让用户决定.
//      另外对于特殊的应用 dna序列 来说,它的k一般比较小,所以4^k也可以提供一个完整的全域向量,至于每一个k_mer的位置编码也不需要
//      通过遍历或者哈希等方式来决定,直接将 dna_shingling.bitset_value cast to size_t 就可以提供一个位置,
//      而且这个位置只要 dna_shingling 的字符串内容确定就一定是唯一的. 如果你用哈希或者遍历方式确定位置编码的话,
//      你换一个哈希算法或者换一种遍历次序得到的位置就可能不同了.
//      因此我还在思考 weight minhash 的第二种实现是否需要改进,因为如果换一种遍历方式(即换一下数据集的数据位置得到的结果就不一样了),
//      但是根据google的论文算法描述,应该是直接遍历一遍非0数据即可,对顺序似乎并没有要求,所以还有待改进...
// TODO
//    Jaccard sim 及其对应的 minhash(Jaccard Index 及其对应的 weight minhash)都是从集合持有元素的角度来衡量相似度的,
//    但是它没有衡量一个序列中 k_mer 的分布是否相似. 仅仅衡量相似的k_mer个数是不够的, 还需要从空间上衡量相似度,
//    这就是 k_mer natural vector 要解决的问题. 可以考虑引入一个权重决定两种相似度的比重,然后综合起来得到最后的相似度.

#endif //LSH_CPP_LSH_COSINE_SIMILARITY_H
