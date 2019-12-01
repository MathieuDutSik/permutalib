#ifndef DEFINE_STAB_CHAIN_MAIN_H
#define DEFINE_STAB_CHAIN_MAIN_H


#include "StabChain.h"
#include "stbcrand.h"


namespace permutalib {

  
// The main function
// Right now we do not implement the PCGS algorithm
//
template<typename Telt, typename Tint>
StabChain<Telt> StabChainOp_listgen(std::vector<Telt> const& Lgen, StabChainOptions<Tint> const& options)
{
  int degree = LargestMovedPoint( Lgen );
  std::cerr << "degree=" << degree << "\n";
  if (degree > 100) {
    std::cerr << "SEARCH : Before call to StabChainRandomPermGroup\n";
    Telt TheId(degree);
    return StabChainRandomPermGroup(Lgen, TheId, options);
  }
  std::cerr << "SEARCH : Doing the ordinary Schreier Sims\n";
  int n=options.n;
  StabChain<Telt> S = EmptyStabChain<Telt>(n);
  if (!IsTrivial_ListGen(Lgen)) {
    S->comm->UseCycle=true;
    std::cerr << "Before call to StabChainStrong\n";
    StabChainStrong(S, Lgen, options );
  }
  if (!options.reduced && options.base.size() > 0) {
    ExtendStabChain(S, options.base);
  }
  /*
    The business with StabChainOptions look eminently dangerous and a reliable replacement
    has to be found.
    It is a record of the options chosen for the stabilizer chain that is outside of the
    variable itself!
  if (options.random > 0) {
        if IsBound( StabChainOptions( Parent( G ) ).random )  then
            options.random := Minimum( StabChainOptions( Parent( G ) ).random,
                                      options.random );
        fi;
        StabChainOptions( G ).random := options.random;
	fi;*/  
  return S;
}



template<typename Telt, typename Tint>
std::pair<bool, StabChain<Telt>> StabChainOp_stabchain(StabChain<Telt> const& G, StabChainOptions<Tint> const& options)
{
  std::cerr << "CPP Before call of StructuralCopy\n";
  StabChain<Telt> S = StructuralCopy(G);
  std::cerr << "CPP After call of StructuralCopy\n";
  if (options.base.size() > 0) {
    if (!ChangeStabChain(S, options.base, options.reduced)) {
      return {false, {}};
    }
  }
  else {
    if (options.reduced) {
      ReduceStabChain(S);
    }
  }
  return {true, S};
}


template<typename Telt, typename Tint>
StabChain<Telt> StabChainOp_stabchain_nofalse(StabChain<Telt> const& G, StabChainOptions<Tint> const& options)
{
  std::cerr << "CPP Before call to StabChainOp_stabchain\n";
  std::pair<bool, StabChain<Telt>> eRec = StabChainOp_stabchain(G, options);
  std::cerr << "CPP After call to StabChainOp_stabchain\n";
  if (!eRec.first) {
    std::cerr << "The nofalse has not been matched\n";
    throw TerminalException{1};
  }
  return eRec.second;
}
 

 

 
template<typename Telt, typename Tint>
Tint Order(StabChain<Telt> const& G)
{
  return SizeStabChain<Telt,Tint>(G);
}

 

template<typename Telt, typename Tint>
StabChain<Telt> MinimalStabChain(std::vector<Telt> const& LGen, int const& n)
{
  StabChainOptions<Tint> options = GetStandardOptions<Tint>(n);
  int largMov=LargestMovedPoint(LGen);
  options.base = ClosedInterval(0, largMov);
  std::cerr << "Before StabChainOp_listgen\n";
  return StabChainOp_listgen(LGen, options);
}





}


#endif
