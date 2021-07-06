#ifndef INCLUDE_GROUP_SERIALIZATION_H
#define INCLUDE_GROUP_SERIALIZATION_H

// Boost serialization

#include <boost/archive/tmpdir.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>




namespace boost::serialization {

  template<class Archive, typename Tgroup>
  inline void load(Archive & ar, Tgroup & val, const unsigned int version) {
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    std::cerr << "load(Tgroup), step 1\n";
    Tidx n_act;
    size_t n_gen;
    ar & make_nvp("n_act", n_act);
    ar & make_nvp("n_gen", n_gen);
    std::vector<Telt> LGen;
    for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
      std::vector<Tidx> eList(n_act);
      for (Tidx i=0; i<n_act; i++) {
        ar & make_nvp("val", eList[i]);
      }
      Telt eGen(eList);
      LGen.push_back(eGen);
    }
    val = Tgroup(LGen, n_act);
    std::cerr << "load(Tgroup), step 2\n";
  }

  template<class Archive, typename Tgroup>
  inline void save(Archive & ar, Tgroup const& val, const unsigned int version) {
    using Telt = typename Tgroup::Telt;
    using Tidx = typename Telt::Tidx;
    std::cerr << "save(Tgroup), step 1\n";
    Tidx n_act = val->n_act();
    ar & make_nvp("n_act", n_act);
    std::vector<Telt> LGen = val.GeneratorsOfGroup();
    size_t n_gen=LGen.size();
    ar & make_nvp("n_gen", n_gen);
    for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
      const Telt& eGen = LGen[i_gen];
      for (Tidx i=0; i<n_act; i++) {
        Tidx val = eGen.at(i);
        ar & make_nvp("val", val);
      }
    }
    std::cerr << "save(Tgroup), step 2\n";
  }

  template<class Archive, typename Tgroup>
  inline void serialize(Archive & ar, Tgroup & val, const unsigned int version) {
    std::cerr << "split_free(Tgroup), step 1\n";
    split_free(ar, val, version);
    std::cerr << "split_free(Tgroup), step 2\n";
  }
}

#endif
