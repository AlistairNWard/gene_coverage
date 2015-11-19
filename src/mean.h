// ***************************************************************************
// Alistair Ward
// Marth Lab, USTAR Center for Genetic Discovery
// University of Utah School of Medicine
// ---------------------------------------------------------------------------
// Last modified: 12 September 2015
// ---------------------------------------------------------------------------
// Calculate mean, median and standard deviation
// ***************************************************************************

#ifndef RUNNING_MEAN_H
#define RUNNING_MEAN_H

#include <vector>

// Data structure for mean calculation.
struct meanData {

  // Set of variables that are reset for each new region.
  double mean;
  double sd;
  double var;
  double M_Old;
  double M_New;
  double S_Old;
  double S_New;
  int k;
  bool isFirst;

  // Set of variables that are not reset and allow calculation of
  // mean across several regions.
  double t_mean;
  double t_sd;
  double t_var;
  double t_M_Old;
  double t_M_New;
  double t_S_Old;
  double t_S_New;
  int t_k;
  bool t_isFirst;

  // Constructor.
  meanData()
    : k(0),
    M_Old(0),
    M_New(0),
    S_Old(0),
    S_New(0),
    isFirst(true),
    mean(0),
    sd(0),
    var(0),
    t_mean(0),
    t_sd(0),
    t_var(0),
    t_M_Old(0),
    t_M_New(0),
    t_S_Old(0),
    t_S_New(0),
    t_k(0),
    t_isFirst(0)
  {}
};

class runningMean {
  public:
    runningMean(void);
    ~runningMean(void);

  // Public methods.
  public:
    struct meanData data;
    void update(int);
    void calculateMean(void);
    void reset(void);
    void resetMultiple(void);
};

#endif // RUNNING_MEAN_H
