#include "hnode.h"

#include "htree.h"
#include "sym_ph.h"
#include "sym_bdexp.h"
#include "sym_state_space_manager.h"
#include "sym_engine.h" 

#include "../debug.h"

#include "original_state_space.h"

using namespace std;

namespace symbolic {

    HNode::HNode(HTree * tree_, const SymParamsMgr & params) : 
	tree(tree_), ph(nullptr),  
	state_space(new OriginalStateSpace(tree->get_engine()->getVars(),
					       params, OperatorCost::NORMAL)){
    }

    HNode::
    HNode(HNode * o, SymPH * ph_, 
	  unique_ptr<SymStateSpaceManager> state_space_mgr) : 
	tree(o->tree), ph(ph_), state_space(std::move(state_space_mgr)){
    }

    void HNode::failed_exploration(SymBDExp * newExp){
	failedForExps.insert(newExp);
	for(auto & p : parents){
	    p->failed_exploration(newExp);
	}
    }

    void HNode::notuseful_exploration(SymBDExp * newExp){
	if(newExp){
	    notUsefulForExps.insert(newExp);
	    notuseful_exploration(newExp->getParent()); //This is not useful for newExp or its parent
	    //All my childs are also not useful
	    for(auto & c : children){
		c->notuseful_exploration(newExp);
	    }
	}
    }

    void HNode::add_exploration(unique_ptr<SymBDExp> && newExp){
	//  if(res->init(this, searchDir, parent, maxTime, maxNodes)){
	exp = move(newExp);
	//parent->setHeuristic(*res, true);
	// cout << "I relaxed the exploration!!" << endl;
	// return true;
	//}else{
	//cout << "I cannot relax the exploration!!!!" << endl;
	//Ensure that we do not add another exploration in this HNode for the same exp
	//failedForExps.insert(parent);
	//return false;
	//}
    }

    bool HNode::hasExpFor(SymBDExp * bdExp) const{
	return (exp && exp->isExpFor(bdExp)) || 
	    failedForExps.count(bdExp) ||
	    notUsefulForExps.count(bdExp);
    }

    bool HNode::isUsefulFor(SymBDExp * bdExp) const {
	return notUsefulForExps.count(bdExp) == 0;
    }

    void HNode::getAllParents(set<HNode *> & setParents){
	for(auto p : parents){
	    setParents.insert(p);
	    p->getAllParents(setParents);
	}
    }

    void HNode::addChildren(HNode * newNode){
	DEBUG_MSG(cout << "ADD CHILDREN IN " << *this << ": " << *newNode << endl;);
	//Create and add the new explorations
	children.push_back(newNode);
    }


    void HNode::addParent(HNode * n){
	//Create and add the new explorations
	parents.push_back(n);
    }


    std::ostream & operator<<(std::ostream &os, const HNode & n){
	os << *n.state_space;

/*os << " with exploration: " <<  *(exp.get()); */
	return os;
    }


    // SymBDExp * HNode::relax(SymBDExp * _exp) const {
    // 	if(ph){
    // 	    return ph->relax(_exp);
    // 	}else{
    // 	    return tree->relax(_exp);
    // 	}
    // }

    SymBDExp * HNode::getPerimeter() const {
	if(!expPerimeters.empty()){  //Use the other perimeter instead
	    DEBUG_PHPDBS(cout << ">> Reusing stored perimeter"<< endl;);
	    return expPerimeters[0].get(); 
	}

	return exp.get();
    }

    SymBDExp * HNode::initSearch(const SymParamsSearch & searchParams, Dir dir) {
	assert(!exp);
	if (!exp) {
	    exp.reset(new SymBDExp (tree->get_engine(), searchParams, dir));
	
	    if(exp->initFrontier(state_space.get(), numeric_limits<int>::max(), numeric_limits<int>::max()) &&
	       exp->initAll(numeric_limits<int>::max(), numeric_limits<int>::max())){
	    }else{
		cout << "Init exploration failed" << endl;
		exit(-1);
	    }
	}
	return exp.get();
    }
}
