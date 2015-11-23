// ***************************************************************************
// Alistair Ward
// Marth Lab, USTAR Center for Genetic Discovery
// University of Utah School of Medicine
// ---------------------------------------------------------------------------
// Last modified: 12 September 2015
// ---------------------------------------------------------------------------
// Calculate mean, median and standard deviation
// ***************************************************************************

#ifndef DATA_PROCESSING_H
#define DATA_PROCESSING_H

#include <vector>

class coverageData {

  public:
    // Hold the data for generating the data.
    std::vector<int> depthData;
    std::vector<int> featureLengths;
    std::vector<int> geneDepth;
  
    // Set of variables that are reset for each new region.
    std::vector<double> featureMean;
    std::vector<double> featureMedian;
    std::vector<double> featureSd;
    std::vector<int> featureMin;
    std::vector<int> featureMax;

    // And the same values for the gene level.
    double geneMean;
    double geneMedian;
    double geneSd;
    int geneMin;
    int geneMax;

  public:
    coverageData(std::size_t);
    ~coverageData(void);

  // Public methods.
  public:
    void update(int);
    void sort();
    float unsortedMedian();
};

#endif // DATA_PROCESSING_H
