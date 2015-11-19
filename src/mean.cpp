// ***************************************************************************
// Alistair Ward
// Marth Lab, USTAR Center for Genetic Discovery
// University of Utah School of Medicine
// ---------------------------------------------------------------------------
// Last modified: 12 September 2015
// ---------------------------------------------------------------------------
// Calculates mean, median and standard deviation
// ***************************************************************************

#include "math.h"
#include "mean.h"
#include <iostream>
using namespace std;

// Constructor
runningMean::runningMean(void) {
}

runningMean::~runningMean(void) {
}

void runningMean::update(int value) {
  
  // Increment the values.
  data.M_Old   = data.M_New;
  data.S_Old   = data.S_New;
  data.t_M_Old = data.t_M_New;
  data.t_S_Old = data.t_S_New;

  // If this is the first value, initialise.
  if (data.t_isFirst == true) {
    data.t_isFirst = false;
    data.t_M_New   = value;
    data.t_S_New   = 0;
    data.t_k       = 1;

  //Otherwise, use the recursion relations.
  } else {
    data.t_k++;
    data.t_M_New = data.t_M_Old + (value - data.t_M_Old) / data.t_k;
    data.t_S_New = data.t_S_Old + (value - data.t_M_Old) * (value - data.t_M_New);
  }

  // If this is the first value, initialise.
  if (data.isFirst == true) {
    data.isFirst = false;
    data.M_New   = value;
    data.S_New   = 0;
    data.k       = 1;

  //Otherwise, use the recursion relations.
  } else {
    data.k++;
    data.M_New = data.M_Old + (value - data.M_Old) / data.k;
    data.S_New = data.S_Old + (value - data.M_Old) * (value - data.M_New);
  }
}

// Calculate the mean and standard deviation.
void runningMean::calculateMean(void) {

  // Calculate the mean and standard deviation for this region.
  data.mean = data.M_New;
  if (data.k > 1) { data.var = data.S_New / (data.k - 1); }
  else { data.var = data.S_Old; }
  data.sd = sqrt(data.var);

  // And for the sum of all regions.
  data.t_mean = data.t_M_New;
  if (data.t_k > 1) { data.t_var = data.t_S_New / (data.t_k - 1); }
  else { data.t_var = data.t_S_Old; }
  data.t_sd = sqrt(data.t_var);
}

// Reset the variables for the individual region mean.
void runningMean::reset(void) {
  data.mean = 0;
  data.sd = 0;
  data.var = 0;
  data.M_Old = 0;
  data.M_New = 0;
  data.S_Old = 0;
  data.S_New = 0;
  data.k = 0;
  data.isFirst = true;
}

// Reset the variables for the mean of multiple regions.
void runningMean::resetMultiple(void) {
  data.t_mean = 0;
  data.t_sd = 0;
  data.t_var = 0;
  data.t_M_Old = 0;
  data.t_M_New = 0;
  data.t_S_Old = 0;
  data.t_S_New = 0;
  data.t_k = 0;
  data.t_isFirst = true;
}
