#include "SortAlgo.h"
#include "Utils.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

// Two-moment bracketing + binary search sort with statistical median and variance
void TwoMomentSort::sort(SortGame &game) {
    // clear previous results
    indexes.clear();
    output_rank.clear();

    int n = game.getSize();
    // Unify predictions and generate statistical predictions (median and variance)
    new_pred();
    generate_statistical_preds();
    
    // Map from predicted ranks to actual indices
    std::vector<int> p_to_A(n);
    for (int i = 0; i < n; ++i) {
        p_to_A[uni_preds[i]] = i;
    }
    
    // Median estimates are already computed in generate_statistical_preds
    
    // Compute bracket widths based on statistical standard deviation
    std::vector<double> sigma(n);
    for (int r = 0; r < n; ++r) {
        // Use the variance computed for the element that will be inserted at position r
        int idx = p_to_A[r];
        sigma[r] = std::sqrt(var_preds[idx]);
    }

    // finger tree setup
    ScapegoatTree T;
    Node* finger = nullptr;
    Comp clean = [&](int a, int b) { return game.compare(a, b); };
    Comp dirty = [&](int a, int b) {
        // inverse compare if involving finger
        if (finger && (a == finger->value || b == finger->value))
            return !clean(a, b);
        return clean(a, b);
    };

    // insert elements in order of predicted rank
    for (int r = 0; r < n; ++r) {
        int ai = p_to_A[r];
        // first element
        if (r == 0) {
            finger = T.root = new Node(ai);
            continue;
        }
        
        // Use our statistical median estimate for this element
        // Add 1.0 because insertion positions are 1-indexed in the tree
        double m = median_estimates[ai] + 1.0;
        double v = sigma[r];
        int Tsz = T.getSize(T.root);

        int l, h;
        // Unified exponential bracket expansion (handles var_preds == 0 as well)
        int d = 1;
        while (true) {
            l = std::max(1, (int)std::floor(m - d * v));
            h = std::min(Tsz, (int)std::ceil(m + d * v));
            // Check if element at l-1 is smaller than ai (proper left bracket)
            bool leftOk  = (l == 1)  || game.compare(T.find_kth_smallest(T.root, l - 1)->value, ai);
            // Check that element at h is not earlier than ai (ai <= element_h)
            bool rightOk = (h == Tsz) || !game.compare(T.find_kth_smallest(T.root, h)->value, ai);
            if (leftOk && rightOk)
                break;                  // bracket correctly encloses ai
            if (l == 1 && h == Tsz)
                break;                  // cannot expand further
            d <<= 1;                    // double the bracket size
        }
        
        // Lower-bound binary search between [l, h]
        int lo = l, hi = h + 1;           // hi is exclusive upper bound (h+1 serves as sentinel)
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (mid == h + 1) {
                lo = mid;
                break;
            }
            Node* midNode = T.find_kth_smallest(T.root, mid);
            if (game.compare(midNode->value, ai))
                lo = mid + 1;           // ai is after midNode
            else
                hi = mid;               // ai is before or equal midNode
        }
        int pos = lo;  // insertion index (1-indexed)
        
        // Rank-based insertion using split/merge – avoids extra comparisons
        {
            auto parts = T.split_root(pos - 1);          // first pos-1 nodes to the left
            Node* mid = new Node(ai);
            T.root = T.join(parts.first, mid, parts.second);
            finger = mid;   // optional: keep for symmetry, not used in compares
        }
    }

    // output sorted order
    T.LinearOutput(&indexes);
    index_to_rank();
}
