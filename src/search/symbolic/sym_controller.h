#ifndef SYMBOLIC_CONTROLLER_H
#define SYMBOLIC_CONTROLLER_H

//Shared class for SymEngine and smas_heuristic

#include "sym_state_space_manager.h"
#include "sym_enums.h"
#include "sym_params_search.h"

#include <vector>
#include <memory>
#include <limits>

namespace options { 
    class OptionParser;
    class Options;
}

using options::OptionParser;
using options::Options;

namespace symbolic {
class SymAbstraction;
class SymSolution;
class SymVariables;
class SymPH;
class SymBDExp;


class SymController { 
 protected:
  SymVariables * vars; //The symbolic variables are declared here  

  SymParamsMgr mgrParams; //Parameters for SymStateSpaceManager configuration.
  SymParamsSearch searchParams; //Parameters to search the original state space 
  
 public:
  SymController(const Options &opts);
  virtual ~SymController() {}

  virtual void new_solution(const SymSolution & /*sol*/){}
  virtual void setLowerBound(int /*lower*/){}
  virtual int getUpperBound() const{return std::numeric_limits<int>::max();}
  virtual int getLowerBound() const{return 0;}
  virtual bool solved () const{return false;}
  virtual SymBDExp * relax(SymBDExp * /*exp*/) const {return nullptr;}

  //HNode * createHNode(HNode * node, std::unique_ptr <SymAbstraction> && abs, SymPH * ph);

  inline SymVariables * getVars(){
    return vars;
  }

  inline const SymParamsMgr & getMgrParams() const{
    return mgrParams;
  }

  inline const SymParamsSearch & getSearchParams() const{
    return searchParams;
  }

  static void add_options_to_parser(OptionParser &parser,
				    int maxStepTime, int maxStepNodes);
};
}
#endif
