
# Distributional Predictions for Sorting
This repository contains our experiments for our CS 2241 class project on distributional predictions for sorting. This codebase is a fork of the official implementation of *Sorting with Predictions*. (Thank you!)

# Run the code
To run the code, please run the following commands in the terminal:
```
bash compile.sh
./main
```
Then, input the name and parameters of the settings you want to test. For example,
```
positional class 10000 10 # positional predictions, class setting, with 10000 items, repeated 10 times.
positional decay 1000 30 # positional predictions, decay setting, with 1000 items, repeated 30 times.
dirty good-dominating 10000 10 # dirty comparison predictions, good-dominating setting, with 10000 items, repeated 10 times.
```

To visualize the number of comparisons needed by each algorithm and baseline, check the Python script `plot.ipynb`.

# Proposed algorithms
- The implementation of Two Moment Sort is at `algorithms/TMS.cpp`.
- The implementation of Distributional Bisection Sort is at `algorithms/DDS.cpp`.
