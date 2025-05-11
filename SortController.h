#pragma once
#include "settings/Interface.h"
#include "algorithms/SortAlgo.h"
#include <vector>

using namespace std;

// External variables declared in main.cpp
extern vector<double> var_preds;
class SortController {
public:
    SortController(SortGame* game);
    void addAlgorithm(SortAlgorithm* algorithm, string name);
    vector<ll> runGame(bool verbose = false);

private:
    SortGame* game;
    vector<SortAlgorithm *> algorithms;
    vector<string> names;
    
    double data_time;
    vector<double> algo_time;
};
