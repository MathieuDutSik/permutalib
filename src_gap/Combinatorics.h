#ifndef DEFINE_PERMUTALIB_COMBINATORICS_H
#define DEFINE_PERMUTALIB_COMBINATORICS_H







template<typename T1, typename T2>
void SortParallel_PairList(std::vector<T1> & list1, std::vector<T2> & list2)
{
  struct PairData {
    T1 val1;
    T2 val2;
  };
  int len=list1.size();
  std::vector<PairData> ListPair(len);
  for (int i=0; i<len; i++)
    ListPair[i] = {list1[i], list2[i]};
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const & a, PairData const& b) -> bool {
         return a.val1 < b.val1;
       });
  for (int i=0; i<len; i++) {
    list1[i] = ListPair[i].val1;
    list2[i] = ListPair[i].val2;
  }
}





#endif
