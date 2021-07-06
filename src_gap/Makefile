#PROGRAM	= GapRepresentativeActionOnSet GapStabilizerOnSet TestStabilizerOnSet TestRepresentativeActionOnSet GapEquivalenceOnSet GapCanonicalImage BenchmarkingPermutalib Factorize TestPermutationProduct
#PROGRAM	= EqSymmGroup
#PROGRAM	= TestStabilizerOnSet
#PROGRAM	= GapStabilizerOnSet
#PROGRAM	= GapRepresentativeActionOnSet
#PROGRAM	= TestRepresentativeActionOnSet
#PROGRAM	= TestSymmetricGroup
#PROGRAM	= TestPermutationProduct
#PROGRAM	= GapCanonicalImage
#PROGRAM	= BenchmarkingPermutalib
#PROGRAM	= BenchmarkingPermutalib GapCanonicalImage
#PROGRAM	= Factorize
PROGRAM	= TestSerialization

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++11
#CC		= clang++ -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++11
#CC		= g++ -I $(BOOST_INCDIR)    -std=c++11 -Wall -g -O2
#CC		= g++ -std=c++14 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fstack-protector-all
#CC		= g++ -std=c++14 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fsanitize=bool -fstack-protector-all
#CC		= g++ -std=c++17 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fstack-protector-all
#CC		= g++ -std=c++17 -Wall -g -O2 -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment -Wstrict-aliasing
#CC		= g++ -std=c++17 -Wall -g -O3
#CC		= g++ -std=c++17 -Wall -Wconversion -g -O3
CC		= g++ -std=c++17 -Wall -Wextra -g -O3
#CC		= clang++ -std=c++17 -Wall -Wconversion -g -O3
#CC		= clang++ -std=c++17 -O3 -fprofile-instr-generate
#CC		= clang++ -std=c++17 -O3 -fprofile-instr-use=a.profdata
#CC		= clang++ -std=c++17 -O3

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -DPERMUTALIB_BLOCKING_SANITY_CHECK
#CLOCAL          = -DPERMUTALIB_BLOCKING_SANITY_CHECK -DSYNCHRONIZED_DEBUG_GAP478

#CFLAGS 		= -I$(BOOST_INCDIR) -I$(GMP_INCDIR) 
CFLAGS 		= -I$(BOOST_INCDIR)

LIBS		= $(GMP_CXX_LINK)

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
