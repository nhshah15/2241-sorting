#include "SortAlgo.h"
#include "Utils.h"
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

// Simplified Two-moment bracketing + binary search sort
void TwoMomentSort_small::sort(SortGame &game) {
    // clear previous results
    indexes.clear();
    output_rank.clear();

    int n = game.getSize();
    // Unify predictions and generate statistical predictions (median and variance)
    new_pred();
    generate_statistical_preds();
    
    // Map from predicted ranks to actual indices
    vector<int> p_to_A(n);
    for (int i = 0; i < n; ++i) {
        p_to_A[uni_preds[i]] = i;
    }
    
    // Median estimates are already computed in generate_statistical_preds
    
    // Compute bracket widths based on statistical standard deviation
    vector<double> sigma(n);
    for (int r = 0; r < n; ++r) {
        // Use the variance computed for the element that will be inserted at position r
        int idx = p_to_A[r];
        sigma[r] = sqrt(var_preds[idx]);
    }

    // Instead of using a tree, we'll use a simple vector to store sorted elements
    vector<int> sorted_elements;
    
    // Insert elements in order of predicted median rank
    for (int r = 0; r < n; ++r) {
        int ai = p_to_A[r];
        
        // First element, just add it
        if (r == 0) {
            sorted_elements.push_back(ai);
            continue;
        }

        // Use our statistical median estimate for this element
        // Add 1.0 because insertion positions are 0-indexed in the vector
        double m = median_estimates[ai] + 1.0;
        double v = sigma[r];
        int sz = sorted_elements.size();
        
        // Determine bracket boundaries, similar to TwoMomentSort
        int l, h;
        int d = 1;
        
        while (true) {
            l = max(0, (int)floor(m - d * v));
            h = min(sz - 1, (int)ceil(m + d * v));
            
            // Check brackets, converted to use array indexes instead of tree nodes
            bool leftOk = (l == 0) || game.compare(sorted_elements[l-1], ai);
            bool rightOk = (h == sz - 1) || !game.compare(sorted_elements[h], ai);
            
            if (leftOk && rightOk)
                break;                  // bracket correctly encloses ai
            if (l == 0 && h == sz - 1)
                break;                  // cannot expand further
            d <<= 1;                    // double the bracket size
        }
        
        // Binary search to find insertion position
        int lo = l, hi = h + 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (mid == sz) {
                lo = mid;
                break;
            }
            if (game.compare(sorted_elements[mid], ai))
                lo = mid + 1;           // ai is after element at mid
            else
                hi = mid;               // ai is before or equal to element at mid
        }
        
        // Insert at the correct position
        sorted_elements.insert(sorted_elements.begin() + lo, ai);
    }
    
    // Copy the final sorted order to indexes
    indexes = sorted_elements;
    index_to_rank();
}
