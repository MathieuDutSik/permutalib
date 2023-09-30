// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Permutation.h"
#include "gmpxx.h"
#include <fstream>

#include "Group.h"

int main(int argc, char *argv[]) {
  try {
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    if (argc != 2) {
      std::cerr << "TestSerialization [EXMP]\n";
      throw permutalib::PermutalibException{1};
    }
    std::string InputFile = argv[1];
    //
    std::ifstream is(InputFile);
    int nbGroup;
    is >> nbGroup;
    for (int iGroup = 0; iGroup < nbGroup; iGroup++) {
      std::cerr << "iGroup=" << iGroup << " / " << nbGroup << "\n";
      Tgroup eG = permutalib::ReadGroupFromStream<Tgroup>(is);
      //
      // Setting the archives
      //
      std::string filename =
          "/tmp/Group_boost_archive_" + std::to_string(iGroup);
      // save data to archive
      {
        std::ofstream ofs(filename);
        boost::archive::text_oarchive oa(ofs);
        // write class instance to archive
        oa << eG;
        // archive and stream closed when destructors are called
      }
      // load data from archive
      Tgroup fG;
      {
        // create and open an archive for input
        std::ifstream ifs(filename);
        boost::archive::text_iarchive ia(ifs);
        // read class state from archive
        ia >> fG;
        // archive and stream closed when destructors are called
      }
      bool test = eG == fG;
      if (!test) {
        std::cerr << "The serialization failed\n";
        throw permutalib::PermutalibException{1};
      }
    }
    std::cerr << "Correct termination of the program\n";
  } catch (permutalib::PermutalibException const &e) {
    std::cerr << "Erroneous completion of the program\n";
    exit(e.eVal);
  }
}
