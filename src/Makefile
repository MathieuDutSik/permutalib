PROGRAM	= GapRepresentativeActionOnSet GapStabilizerOnSet GapEquivalenceOnSet GapCanonicalImage BenchmarkingPermutalib Factorize TestPermutation TestSymmetricGroup GroupProperties EqSymmGroup TestSerialization GapDerivedSubgroup GapSmallGeneratingSet TestEnumerationElement GapCentreSubgroup GapCentralizerElt GapComputeSequenceBlockSystems GapAscendingChain GapAscendingChainSubgroup TestRightCosets TestLeftCosets TestComputeDoubleCosets GapComputeDoubleCosets TestComputeDoubleCosetsStabilizers TestPreImageSubgroup
#PROGRAM	= EqSymmGroup
#PROGRAM	= GapStabilizerOnSet
#PROGRAM	= GapRepresentativeActionOnSet
#PROGRAM	= TestSymmetricGroup
#PROGRAM	= TestPermutation BenchmarkingPermutalib
#PROGRAM	= GapCanonicalImage
#PROGRAM	= GroupProperties
#PROGRAM	= GapSmallGeneratingSet
#PROGRAM	= GapDerivedSubgroup
#PROGRAM	= GapCentreSubgroup
#PROGRAM	= GapCentralizerElt
#PROGRAM	= Factorize
#PROGRAM	= TestSerialization
#PROGRAM	= TestEnumerationElement
#PROGRAM	= GapComputeSequenceBlockSystems
#PROGRAM	= GapAscendingChain
#PROGRAM	= GapIntersection
#PROGRAM	= BenchmarkingPermutalib

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++17 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++17
#CC		= clang++ -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++17
#CC		= g++ -I $(BOOST_INCDIR)    -std=c++17 -Wall -g -O2
#CC		= g++ -std=c++17 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fsanitize=bool -fstack-protector-all
#CC		= g++ -std=c++17 -Wall -g -O3 -fsanitize=address -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fstack-protector-all
CC		= g++ -std=c++17 -Wall -Wextra -O3
#CC		= g++ -std=c++17 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fstack-protector-all
#CC		= g++ -std=c++17 -Wall -g -O2 -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment -fsanitize=leak -Wstrict-aliasing
#CC		= g++ -std=c++17 -Wall -Wextra -O3
#CC		= g++ -pg -std=c++17 -Wall -g -O3
#CC		= g++ -std=c++17 -Wall -g -O3
#CC		= g++ -std=c++17 -Wall -Wconversion -g -O3
#CC		= g++ -std=c++17 -Wall -Wextra -g -O3
#CC		= clang++ -std=c++17 -Wall -Wconversion -g -O3
#CC		= clang++ -std=c++17 -O3 -fprofile-instr-generate
#CC		= clang++ -std=c++17 -O3 -fprofile-instr-use=a.profdata
#CC		= clang++ -std=c++17 -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

#CLOCAL          =
CLOCAL          = -DPERMUTALIB_BLOCKING_SANITY_CHECK -DCHECK_COSET_ENUMERATION -DDEBUG_PRE_IMAGE_COMPLEXITY -DTIMINGS_PRE_IMAGE -DTIMINGS_SMALL_GENERATING_SET -DDEBUG_SMALL_GENERATING_SET
#CLOCAL          = -DPERMUTALIB_BLOCKING_SANITY_CHECK -DSYNCHRONIZED_DEBUG_GAP478
#CLOCAL          = -DPERMUTALIB_BLOCKING_SANITY_CHECK -DSYNCHRONIZED_DEBUG_GAP478 -DDEBUG -DDEBUG_ADD_GEN_SCH -DDEBUG_CHANGE_STAB_CHAIN -DDEBUG_EQUALITY -DDEBUG_INV_REP -DDEBUG_STABCHAIN







CFLAGS 		= -I$(BOOST_INCDIR) -I$(GMP_INCDIR) 

LIBS		= $(GMP_CXX_LINK) $(BOOST_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CFLAGS) $(CLOCAL) $(NOLINKFLAG) $<



#	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
