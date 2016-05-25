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

#include <string>
#include <sstream>
#include <vector>

using namespace std;

class coverageData {

  public:

    // id for the region.
    std::vector<string> ids;

    // Hold the data for generating the data.
    std::vector<int> featureLengths;
  
    // Set of variables that are reset for each new region.
    std::vector<double> featureMean;
    std::vector<double> featureMedian;
    std::vector<double> featureQ1;
    std::vector<double> featureQ3;
    std::vector<double> featureIqr;
    std::vector<double> featureSd;
    std::vector<int> featureMin;
    std::vector<int> featureMax;

    // And the same values for the gene level.
    vector<int> geneCoverage;
    double geneMean;
    double geneMedian;
    double geneQ1;
    double geneQ3;
    double geneIqr;
    double geneSd;
    int geneMin;
    int geneMax;

  public:
    coverageData(std::size_t);
    ~coverageData(void);

  // Public methods.
  public:
    void noCoverage();
    void processFeature(vector<int>&, int);
    void processGene();
};

#endif // DATA_PROCESSING_H
