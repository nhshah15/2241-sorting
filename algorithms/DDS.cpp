#include "SortAlgo.h"
#include "Utils.h"
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <memory>
#include <map>
using namespace std;

// Global variables to track DDS path usage across runs
int dds_concentrated_total = 0;
int dds_entropy_guided_total = 0;
int dds_total_elements = 0;

// Extern declaration of prediction-related vectors
extern vector<double> var_preds;
extern vector<double> median_estimates;

/**
 * Structure to hold data for learning-based PMF generation
 */
struct ObservationData {
    // Relative position where elements were actually placed
    vector<double> observed_relative_positions;
    // Prediction errors (absolute difference between predicted and actual position)
    vector<double> prediction_errors;
    // Map of element characteristics to actual positions
    map<int, vector<int>> characteristic_to_positions;
    
    void addObservation(double predicted_pos, double actual_pos, int characteristics) {
        // Store relative position (actual_pos - predicted_pos) / range
        double relative_pos = actual_pos - predicted_pos;
        observed_relative_positions.push_back(relative_pos);
        
        // Store absolute prediction error
        prediction_errors.push_back(abs(relative_pos));
        
        // Store position by element characteristics (could be based on element value, etc.)
        characteristic_to_positions[characteristics].push_back(actual_pos);
    }
    
    // Simple weighted KDE for probability estimation
    double estimateProbability(double x, double bandwidth) {
        if (observed_relative_positions.empty()) return 0.0;
        
        double prob = 0.0;
        for (double obs : observed_relative_positions) {
            // Gaussian kernel
            double z = (x - obs) / bandwidth;
            prob += exp(-0.5 * z * z);
        }
        
        return prob / observed_relative_positions.size();
    }
};

// Global observation data for learning from previous insertions
ObservationData dds_observation_data;

/**
 * Determines if the predicted position should be treated as concentrated
 * @param variance Variance of the prediction
 * @return True if the variance is small enough to be considered concentrated
 */
bool isConcentratedPrediction(double variance) {
    return variance < 0.5; // Simple threshold for high confidence
}

/**
 * Finds the concentrated position if this is a concentrated distribution
 * @param pmf The probability mass function
 * @return Position index if concentrated, -1 otherwise
 */
int findConcentratedPosition(const vector<double>& pmf) {
    auto it = find(pmf.begin(), pmf.end(), 1.0);
    if (it != pmf.end()) {
        return (int)(it - pmf.begin()) + 1;  // Convert to 1-based position
    }
    return -1;  // Not a concentrated distribution
}

/**
 * Performs adaptive bracketing to find search range
 * @param game The sort game instance
 * @param T The scapegoat tree
 * @param ai Current element index
 * @param m Predicted position
 * @param sigma Standard deviation of prediction
 * @return Pair of lo and hi bounds for the search
 */
pair<int, int> adaptiveBracketing(SortGame &game, ScapegoatTree &T, int ai, int m, double sigma) {
    int Tsz = T.getSize(T.root);
    int d = 1;
    int lo = max(1, (int)floor(m - d * sigma));
    int hi = min(Tsz + 1, (int)ceil(m + d * sigma));
    
    // Expand brackets until we confirm they contain the insertion point
    while (true) {
        // Check if element at lo-1 is smaller than ai (proper left bracket)
        bool leftOk = (lo == 1) || game.compare(T.find_kth_smallest(T.root, lo - 1)->value, ai);
        // Check that element at hi is not earlier than ai (ai <= element_hi)
        bool rightOk = (hi == Tsz + 1) || !game.compare(T.find_kth_smallest(T.root, hi)->value, ai);
        
        if (leftOk && rightOk)
            break;  // bracket correctly encloses ai
        if (lo == 1 && hi == Tsz + 1)
            break;  // cannot expand further
            
        d <<= 1;  // double the bracket size
        lo = max(1, (int)floor(m - d * sigma));
        hi = min(Tsz + 1, (int)ceil(m + d * sigma));
    }
    
    return {lo, hi};
}

// Distributional Displacement Sort (DDS) - Strictly following pseudocode
void DDS::sort(SortGame &game) {
    int n = game.getSize();
    
    // Initialize data structures
    indexes.clear();
    output_rank.clear();
    vector<int> sorted_indexes; // T in the pseudocode (the sorted sequence so far)
    
    // Get unique predictions
    new_pred();
    
    // Calculate variance for each element
    var_preds.resize(n);
    for (int i = 0; i < n; ++i) {
        int predicted_rank = uni_preds[i];
        double error_fraction = double(abs(predicted_rank - i)) / n;
        var_preds[i] = error_fraction * error_fraction * n;
        var_preds[i] = std::max(0.5, var_preds[i]);
    }
    
    // Process each element in turn
    for (int i = 0; i < n; ++i) {
        // Pre-insertion check for first element
        if (i == 0) {
            sorted_indexes.push_back(i); // Insert first element at position 1
            continue;
        }
        
        // Get predicted rank and variance for current element
        double pred_rank = uni_preds[i];
        double variance = var_preds[i];
        double sigma = sqrt(variance);
        
        // Generate PMF (simplified for performance)
        bool is_concentrated = (variance < 0.5); // Concentrated distribution check
        int insertion_pos;
        
        if (is_concentrated) {
            // Special case: concentrated distribution
            int m = round(pred_rank + 1.0); // Convert to 1-based position
            m = std::clamp(m, 1, (int)sorted_indexes.size()+1);
            
            // Binary search around position m
            int lo = 0;
            int hi = sorted_indexes.size();
            
            while (lo < hi) {
                int mid = lo + (hi - lo) / 2;
                if (game.compare(i, sorted_indexes[mid])) {
                    // Element i comes before element at mid
                    hi = mid;
                } else {
                    // Element i comes after element at mid
                    lo = mid + 1;
                }
            }
            
            insertion_pos = lo;
        } else {
            // Entropy-guided search on PMF support as per pseudocode
            int ell = 0; // Lower bound (0-based index, equivalent to 1 in pseudocode)
            int h = sorted_indexes.size(); // Upper bound (equivalent to min(i, T.size()) in pseudocode)
            
            // Generate a simplified Gaussian PMF centered at predicted rank
            vector<double> pmf(h - ell + 1, 0.0);
            double center = pred_rank; // 0-based predicted position
            
            // Only calculate PMF values within a reasonable window to save computation
            int window_size = min(h - ell, 5 + (int)(2 * sigma));
            int pmf_start = max(0, (int)round(center) - window_size/2);
            int pmf_end = min(h - ell, (int)round(center) + window_size/2);
            
            // Calculate PMF within window
            for (int k = pmf_start; k <= pmf_end; k++) {
                double z = (k - center) / sigma;
                pmf[k] = exp(-0.5 * z * z);
            }
            
            // Normalize the PMF
            double sum_pmf = 0.0;
            for (auto p : pmf) sum_pmf += p;
            if (sum_pmf > 0) {
                for (auto &p : pmf) p /= sum_pmf;
            }
            
            // Main entropy-guided search loop
            while (ell < h) {
                // Calculate the total probability mass in the current interval
                double total_mass = 0.0;
                for (int k = ell; k < h && k < pmf.size(); k++) {
                    total_mass += pmf[k];
                }
                
                if (total_mass <= 0.0) {
                    // Fallback to binary search when no probability mass
                    int r_star = ell + (h - ell) / 2;
                    
                    if (r_star >= sorted_indexes.size() || game.compare(i, sorted_indexes[r_star])) {
                        h = r_star;
                    } else {
                        ell = r_star + 1;
                    }
                } else {
                    // Find the split point that divides probability mass in half
                    double cumulative = 0.0;
                    int r_star = ell;
                    
                    for (int k = ell; k < h && k < pmf.size(); k++) {
                        cumulative += pmf[k];
                        if (cumulative >= total_mass / 2.0) {
                            r_star = k;
                            break;
                        }
                    }
                    
                    // Ensure r_star is within bounds
                    r_star = max(ell, min(h-1, r_star));
                    
                    // Perform the comparison and narrow search interval
                    if (r_star >= sorted_indexes.size() || game.compare(i, sorted_indexes[r_star])) {
                        // Current element comes before element at r_star
                        h = r_star;
                    } else {
                        // Current element comes after element at r_star
                        ell = r_star + 1;
                    }
                    
                    // Renormalize PMF on the new interval (as per pseudocode)
                    if (ell < h) {
                        double interval_mass = 0.0;
                        for (int k = ell; k < h && k < pmf.size(); k++) {
                            interval_mass += pmf[k];
                        }
                        
                        if (interval_mass > 0.0) {
                            for (int k = ell; k < h && k < pmf.size(); k++) {
                                pmf[k] /= interval_mass;
                            }
                        }
                    }
                }
            }
            
            insertion_pos = ell;
        }
        
        // Insert element i at the determined position
        sorted_indexes.insert(sorted_indexes.begin() + insertion_pos, i);
    }
    
    // Build final output
    indexes.assign(sorted_indexes.begin(), sorted_indexes.end());
    
    // Convert to output_rank
    output_rank.resize(n);
    for (int i = 0; i < n; i++) {
        output_rank[indexes[i]] = i;
    }
}
