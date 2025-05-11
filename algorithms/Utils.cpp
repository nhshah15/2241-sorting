#include "Utils.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <random>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
using namespace std;

vector<double> median_estimates;

void new_pred() {
    // input: preds
    // output: uni_preds
    int n = preds.size();

    buckets.resize(n + 1);
    for (int i = 0; i <= n; i++)
        buckets[i].clear();
    uni_preds.resize(n);
    
    for (int i = 0; i < n; i++)
        buckets[preds[i]].push_back(i);
    
    //random shuffle each bucket
    for (int i = 0; i <= n; i++)
        for (int j = 0; j < buckets[i].size(); j++)
            swap(buckets[i][j], buckets[i][rand() % buckets[i].size()]);

    int counter = 0;
    for (int i = 0; i <= n; i++)
        for (int j = 0; j < buckets[i].size(); j++)
            uni_preds[buckets[i][j]] = counter++;
    
}

void index_to_rank() {
    // input: indexes
    // output: output_rank
    output_rank.clear();
    output_rank.resize(indexes.size());
    for (int i = 0; i < indexes.size(); i++) {
        output_rank[indexes[i]] = i;
    }
}

void output_to_file(vector<string> names, vector<vector<vector<ll>>> results, string filename) {
    //output to a file
    int n = results.size();
    int m = results[0].size();
    // cerr << "sizes" << names.size() << " " << n << " " << m << endl;
    for (int i = 0; i < m; i++) {
        //align
        cout << names[i] << " ";
        for (int j = 0; j < n; j++) {
            vector <ll> reps = results[j][i];
            // cout << names[i] << " ";
            ll sum = 0;
            for (int k = 0; k < reps.size(); k++)
            {
                // cout << reps[k] << " ";
                sum += reps[k];
            }
            cout << sum / reps.size() << " ";
        }
        cout << endl;
    }

    ofstream outFile("data/" + filename + ".txt");
    cerr << "saved to " << filename << endl;
    for (int i = 0; i < m; i++) {
        outFile << names[i] << " ";
    }
    outFile << endl
            << m << " " << n << " " << results[0][0].size() << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            vector <ll> reps = results[j][i];
            for (int k = 0; k < reps.size(); k++) {
                outFile << reps[k] << " ";
            }
            outFile << endl;
        }
    }
    outFile.close();
}

double get_time() {
    return clock() / (double)CLOCKS_PER_SEC;
}


typedef long long ll;
pair<vector<ll>, vector<ll>> parsePopulationData(const string& filename, int old) {
    ifstream file(filename);
    string line;
    vector<ll> row2021, rowOld;
    int targetYear = 2021 - old;

    getline(file, line);
    getline(file, line); //countries

    while (getline(file, line))
    {
        stringstream ss(line);
        int year;
        ss >> year;
        // cout << "see the target year: " << year << endl;
        if (year == 2021 || year == targetYear) {
            vector<ll> rowData;
            ll value;
            while (ss >> value) {
                // cerr << "value: " << value << endl;
                rowData.push_back(value);
            }

            if (year == 2021) {
                row2021 = rowData;
            }
            if (year == targetYear) {
                rowOld = rowData;
            }
        }

        if (!row2021.empty() && !rowOld.empty()) {
            break;
        }
    }

    // output row2021 and rowOld
    // cout << "row2021: ";
    // for (int i = 0; i < row2021.size(); i++) {
    //     cout << row2021[i] << " ";
    // }
    // cout << endl;
    // cout  << "rowOld: ";
    // for (int i = 0; i < rowOld.size(); i++) {
    //     cout << rowOld[i] << " ";
    // }
    // cout << endl;


    return make_pair(row2021, rowOld);
}

vector <vector<bool> > parseTennisRelation(const string &filename, int size) 
{
    ifstream file(filename);
    //first line is the number of players
    int n;
    file >> n;
    cout << "n: " << n << endl;
    //the rest of the file is the relation
    string line;
    vector <vector<bool> > relation;
    getline(file, line);
    while (getline(file, line))
    {
        stringstream ss(line);
        vector<bool> row;
        bool value;
        while (ss >> value)
            row.push_back(value);
        relation.push_back(row);
    }
    
    return relation;
}

void generate_statistical_preds() {
    // Input: uni_preds (after new_pred has been called)
    // Output: var_preds (variance predictions) and median_estimates
    
    int n = uni_preds.size();
    
    // Create a copy of unified predictions to compute statistical properties
    std::vector<int> pred_copies(n);
    for (int i = 0; i < n; ++i) {
        pred_copies[i] = uni_preds[i];
    }
    
    // Calculate median of the predicted distribution
    std::vector<int> sorted_preds = pred_copies;
    std::sort(sorted_preds.begin(), sorted_preds.end());
    double median_pred;
    if (n % 2 == 0) {
        median_pred = (sorted_preds[n/2 - 1] + sorted_preds[n/2]) / 2.0;
    } else {
        median_pred = sorted_preds[n/2];
    }
    
    // Generate median estimates for each element
    median_estimates.resize(n);
    for (int i = 0; i < n; ++i) {
        int predicted_rank = uni_preds[i];
        double deviation_from_median = predicted_rank - median_pred;
        median_estimates[i] = median_pred + deviation_from_median;
    }
    
    // Calculate mean of the predicted distribution
    double mean_pred = 0.0;
    for (const int& pred : pred_copies) {
        mean_pred += pred;
    }
    mean_pred /= n;
    
    // Calculate variance of the predicted distribution
    double variance_pred = 0.0;
    for (const int& pred : pred_copies) {
        double diff = pred - mean_pred;
        variance_pred += diff * diff;
    }
    variance_pred /= n;
    double std_dev_pred = std::sqrt(variance_pred);
    
    // Resize and initialize var_preds
    var_preds.resize(n);
    
    // Calculate prediction quality metrics
    double mean_abs_error = 0.0;
    for (int i = 0; i < n; ++i) {
        mean_abs_error += std::abs(pred_copies[i] - i);
    }
    mean_abs_error /= n;
    
    // Prediction quality factor (0.0 = poor predictions, 1.0 = perfect predictions)
    double prediction_quality = std::max(0.0, 1.0 - mean_abs_error / std::max(1.0, n / 8.0));
    
    // Calculate prediction variance for each element
    for (int i = 0; i < n; ++i) {
        int predicted_rank = uni_preds[i];
        
        // Calculate Z-score (how many std deviations from mean)
        double z_score = (predicted_rank - mean_pred) / std::max(1.0, std_dev_pred);
        
        // Calculate deviation from median
        double deviation_from_median = predicted_rank - median_pred;
        
        // Calculate relative distance with respect to the total range of predictions
        double range = sorted_preds[n-1] - sorted_preds[0];
        double normalized_distance = std::abs(deviation_from_median) / std::max(1.0, range / 8.0);
        
        // Apply exponential decay when predictions are good
        double scaling_exponent = 1.0 - 0.5 * prediction_quality;
        double scaling_factor = std::min(3.0, 1.0 + std::pow(normalized_distance, scaling_exponent));
        
        // Base variance decreases quadratically with prediction quality
        double quality_factor = 1.0 - prediction_quality * prediction_quality;
        double base_variance = std_dev_pred * quality_factor;
        
        // Adaptive cap based on prediction quality
        double max_variance = std::min(n / (4.0 + 12.0 * prediction_quality), 250.0 * (1.0 - 0.8 * prediction_quality));
        
        // Final variance with improved quality-sensitive scaling
        var_preds[i] = std::min(max_variance, base_variance * scaling_factor);
    }
}