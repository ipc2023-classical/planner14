#ifndef SYMBOLICBOLIC_HNODE_H
#define SYMBOLICBOLIC_HNODE_H

#include "sym_state_space_manager.h"
#include <memory> 
#include <vector> 
#include <set> 

namespace symbolic {
    class HTree;
    class SymBDExp;
    class SymParamsSearch;

    class SymPH;
    class SymController;
    class SymExploration;
    class SymParamsMgr;

    class HNode {
    private:
	HTree * tree;
	SymPH * ph;

	std::unique_ptr<SymStateSpaceManager> state_space;

	std::unique_ptr<SymBDExp> exp;

	//If should abstract is activated, we store a list with perimeters
	//of the original exploration to initialize our abstractions. 
	std::vector<std::unique_ptr<SymBDExp>> expPerimeters;

	std::vector <HNode *> children; //Nodes more abstracted
	std::vector <HNode *> parents; //Nodes less abstracted

	std::set <SymBDExp *> failedForExps; //Set of exps we failed to abstract
	std::set <SymBDExp *> notUsefulForExps; //Set of exps we are not useful for
  
    public:
	// Constructor for the original state space
	HNode(HTree * tree, const SymParamsMgr & mgr);
	
	// Constructor for abstract state space
	HNode(HNode * o, SymPH * ph, std::unique_ptr<SymStateSpaceManager> abs); 

	HNode(const HNode & o) = delete;
	HNode(HNode &&) = default;
	HNode& operator=(const HNode& ) = delete;
	HNode& operator=(HNode &&) = default;
	~HNode() = default; 

	SymBDExp * initSearch(const SymParamsSearch & searchParams, Dir dir);

	void getAllParents(std::set<HNode *> & setParents);

	HNode * getParent(){
	    return parents[0];
	}

	void add_exploration(std::unique_ptr<SymBDExp> && newExp);
	void failed_exploration(SymBDExp * newExp);
	void notuseful_exploration(SymBDExp * newExp);

	void addChildren(HNode * newNode);
	void addParent(HNode * newNode);

	bool empty() const{
	    return !exp && failedForExps.empty() && notUsefulForExps.empty();
	}
	bool hasExpFor(SymBDExp * bdExp) const;
	bool isUsefulFor(SymBDExp * bdExp) const;

	inline int numVariablesToAbstract() const {
	    return state_space->numVariablesToAbstract();
	}

	inline int numVariablesAbstracted() const {
	    return state_space->numVariablesAbstracted();
	}

	inline bool isAbstracted() const{
	    return state_space->isAbstracted();
	}

	/*inline const std::vector <HNode *> & getChildren(SymPH * of_ph){
	  return children;
	  }*/

	inline std::vector <HNode *> getChildren(SymPH * of_ph){
	    std::vector <HNode *> res;
	    for(auto c : children){
		if(c->ph == of_ph) 
		    res.push_back(c);
	    }
	    return res;
	}

	/* SymBDExp * relax(SymBDExp * _exp) const; */
  
	SymBDExp * getPerimeter () const;

	inline SymStateSpaceManager * getStateSpace() const{
	    return state_space.get();
	}

	inline SymBDExp * getExp() const{
	    return exp.get();
	}

	inline bool hasStoredPerimeters () const {
	    return !expPerimeters.empty();
	}

	void addPerimeter (std::unique_ptr<SymBDExp> & perimeter) {
	    expPerimeters.push_back(std::move(perimeter));
	} 

	friend std::ostream & operator<<(std::ostream &os, const HNode & n);
    };

}
#endif
