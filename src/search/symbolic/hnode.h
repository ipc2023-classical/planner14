#ifndef SYMBOLIC_HNODE_H
#define SYMBOLIC_HNODE_H

#include "sym_state_space_manager.h"
#include <memory>
#include <vector>
#include <set>

namespace symbolic {
class HTree;
class BDAstar;
class SymParamsSearch;

class SymPH;
class SymController;
class SymExploration;
class SymParamsMgr;

class HNode {
private:
    HTree *tree;
    SymPH *ph;

    std::shared_ptr<SymStateSpaceManager> state_space;

    std::unique_ptr<BDAstar> exp;

    //If should abstract is activated, we store a list with perimeters
    //of the original exploration to initialize our abstractions.
    std::vector<std::unique_ptr<BDAstar>> expPerimeters;

    std::vector <HNode *> children;     //Nodes more abstracted
    std::vector <HNode *> parents;     //Nodes less abstracted

    std::set <BDAstar *> failedForExps;     //Set of exps we failed to abstract
    std::set <BDAstar *> notUsefulForExps;     //Set of exps we are not useful for

public:
    // Constructor for the original state space
    HNode(HTree *tree, const SymParamsMgr &mgr);

    // Constructor for abstract state space
    HNode(HNode *o, SymPH *ph, std::unique_ptr<SymStateSpaceManager> abs);

    HNode(const HNode &o) = delete;
    HNode(HNode &&) = default;
    HNode &operator=(const HNode &) = delete;
    HNode &operator=(HNode &&) = default;
    ~HNode() = default;

    BDAstar *initSearch(const SymParamsSearch &searchParams, Dir dir);

    void getAllParents(std::set<HNode *> &setParents);

    HNode *getParent() {
        return parents[0];
    }

    void add_exploration(std::unique_ptr<BDAstar> newExp);
    void failed_exploration(BDAstar *newExp);
    void notuseful_exploration(BDAstar *newExp);

    void addChildren(HNode *newNode);
    void addParent(HNode *newNode);

    bool empty() const {
        return !exp && failedForExps.empty() && notUsefulForExps.empty();
    }
    bool hasExpFor(BDAstar *bdExp) const;
    bool isUsefulFor(BDAstar *bdExp) const;

    inline int numVariablesToAbstract() const {
        return state_space->numVariablesToAbstract();
    }

    inline int numVariablesAbstracted() const {
        return state_space->numVariablesAbstracted();
    }

    inline bool isAbstracted() const {
        return state_space->isAbstracted();
    }

    /*inline const std::vector <HNode *> & getChildren(SymPH * of_ph){
      return children;
      }*/

    inline std::vector <HNode *> getChildren(SymPH *of_ph) {
        std::vector <HNode *> res;
        for (auto c : children) {
            if (c->ph == of_ph)
                res.push_back(c);
        }
        return res;
    }

    /* BDAstar * relax(BDAstar * _exp) const; */

    BDAstar *getPerimeter() const;

    SymStateSpaceManager *getStateSpace() const {
        return state_space.get();
    }

    std::shared_ptr<SymStateSpaceManager> &getStateSpaceRef() {
        return state_space;
    }

    inline BDAstar *getExp() const {
        return exp.get();
    }

    inline bool hasStoredPerimeters() const {
        return !expPerimeters.empty();
    }

    void addPerimeter(std::unique_ptr<BDAstar> perimeter) {
        expPerimeters.push_back(std::move(perimeter));
    }

    friend std::ostream &operator<<(std::ostream &os, const HNode &n);
};
}
#endif
