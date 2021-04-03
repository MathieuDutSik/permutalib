#PROGRAM	= GapRepresentativeActionOnSet GapStabilizerOnSet TestStabilizerOnSet TestRepresentativeActionOnSet GapEquivalenceOnSet GapCanonicalImage
#PROGRAM	= TestStabilizerOnSet
#PROGRAM	= GapStabilizerOnSet
#PROGRAM	= GapRepresentativeActionOnSet
#PROGRAM	= TestRepresentativeActionOnSet
#PROGRAM	= TestSymmetricGroup
PROGRAM	= GapCanonicalImage

OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2
#CC		= clang++ -Weverything -std=c++11
#CC		= clang++ -Weverything -fsanitize=thread -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fsanitize=bounds -fsanitize=return -fsanitize=null -std=c++11
#CC		= g++ -I $(BOOST_INCDIR)    -std=c++11 -Wall -g -O2
#CC		= g++ -std=c++14 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fstack-protector-all
#CC		= g++ -std=c++14 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fsanitize=bool -fstack-protector-all
#CC		= g++ -std=c++17 -Wall -g -O2 -fsanitize=address -fsanitize=undefined -fsanitize=bool -fsanitize=bounds -fstack-protector-all
CC		= g++ -std=c++17 -Wall -g -O2

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -DPERMUTALIB_BLOCKING_SANITY_CHECK

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
