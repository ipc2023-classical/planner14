#include "bidirectional_search.h"

#include "../utils/debug_macros.h"
#include "sym_engine.h"
#include <algorithm>    // std::reverse
#include <memory>    

#include "../global_operator.h"

using namespace std;
using utils::g_timer;

namespace symbolic {

    BidirectionalSearch::BidirectionalSearch(const SymParamsSearch &params, std::unique_ptr<UnidirectionalSearch> _fw,
					     unique_ptr<UnidirectionalSearch> _bw) : SymSearch(params),
											  fw(std::move(_fw)),
											  bw(std::move(_bw)){
    }

   

    UnidirectionalSearch *BidirectionalSearch::selectBestDirection() const {

	bool fwSearchable = fw->isSearchable();
	bool bwSearchable = bw->isSearchable();
	if (fwSearchable && !bwSearchable) {
	    return fw.get();
	} else if (!fwSearchable && bwSearchable) {
	    return bw.get();
	}

	return fw->nextStepNodes() <= bw->nextStepNodes() ? fw.get() : bw.get();
    }


    bool BidirectionalSearch::finished() const {
	return fw->finished() || bw->finished();
    }


    void BidirectionalSearch::statistics() const {
	if (fw)
	    fw->statistics();
	if (bw)
	    bw->statistics();
	cout << endl;
    }


    bool BidirectionalSearch::stepImage(int maxTime, int maxNodes)  {
	return selectBestDirection()->stepImage(maxTime, maxNodes);
    }

}
