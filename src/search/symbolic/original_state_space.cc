#include "original_state_space.h" 

#include "../debug.h"
using namespace std;

namespace symbolic {

    OriginalStateSpace::OriginalStateSpace(SymVariables * v, 
					   const SymParamsMgr & params, 
					   OperatorCost cost_type_) : 
	SymStateSpaceManager(v, params, cost_type_) {

    }



    void OriginalStateSpace::init_initial_state() {    
	initialState = vars->getStateBDD(g_initial_state_data);
    }

    void OriginalStateSpace::init_goal() {
	goal = vars->getPartialStateBDD(g_goal);
    }

    void OriginalStateSpace::init_mutex(const std::vector<MutexGroup> & mutex_groups){
	//If (a) is initialized OR not using mutex OR edeletion does not need mutex
	if(mutexInitialized || p.mutex_type == MutexType::MUTEX_NOT)
	    return; //Skip mutex initialization
 
	if(p.mutex_type == MutexType::MUTEX_EDELETION){
	    SymStateSpaceManager::init_mutex(mutex_groups, true, true);
	}else{
	    SymStateSpaceManager::init_mutex(mutex_groups, true, false);
	}
    }

    void OriginalStateSpace::init_individual_trs(){
	if(!indTRs.empty()) return;

	DEBUG_MSG(cout << "Initialize individual TRs of original state space" << endl;);
	for(size_t i = 0; i < g_operators.size(); i++){
	    const GlobalOperator * op = &(g_operators[i]);
	    // Skip irrelevant operators 
	    /*if (op->is_dead()){ 
	      continue;
	      }*/
	    int cost = get_adjusted_action_cost(*op, cost_type); 
	    DEBUG_MSG(cout << "Creating TR of op " << i << " of cost " << cost << endl;);
	    indTRs[cost].push_back(move(SymTransition(vars, op, cost)));
	    if(p.mutex_type ==MutexType::MUTEX_EDELETION){
		indTRs[cost].back().edeletion(*this);
	    }
	}
    }

    void OriginalStateSpace::init_transitions(){
	if(!transitions.empty()) return; //Already initialized!
	init_individual_trs(); 
	init_transitions_from_individual_trs();
    }

}
