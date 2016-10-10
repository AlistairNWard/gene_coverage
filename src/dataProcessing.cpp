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
#include "dataProcessing.h"
#include <algorithm>
#include <iostream>
using namespace std;

// Constructor
coverageData::coverageData(size_t size) {

  // Initialise arrays.
  featureLengths[size];
}

coverageData::~coverageData(void) {
}

// If a feature has no coverage, add the statistics to the correct fields.
void coverageData::noCoverage() {
  featureMin.push_back(0);
  featureMax.push_back(0);
  featureMean.push_back(0);
  featureMedian.push_back(0);
  featureQ1.push_back(0);
  featureQ3.push_back(0);
  featureIqr.push_back(0);
  featureSd.push_back(0);
}

// Process a single feature.
void coverageData::processFeature(vector<int>& coverage, int start) {

  // Initialise variables.
  int length  = coverage.size() - start;
  int sum     = 0;
  double sd   = 0.;
  vector<int> newCoverage;

  // Loop over the coverage data for the feature.
  vector<int>::iterator iter    = coverage.begin() + start - 1;
  vector<int>::iterator iterEnd = coverage.end();
  for (; iter != iterEnd; ++iter) {

    // keep a sum of the coverage values.
    sum += *iter;

    // Add the value to the newCoverage vector. This will create the vector with the correct
    // number of entries to be sorted for calculating the median.
    newCoverage.push_back(*iter);
    geneCoverage.push_back(*iter);
  }

  // Calculate the mean.
  double mean = sum / double(length);
  featureMean.push_back(mean);

  // Sort the newCoverage vector.
  std::sort(newCoverage.begin(), newCoverage.end());

  // Now calculate the standard deviation for the feature.
  iter    = newCoverage.begin();
  iterEnd = newCoverage.end();
  for (; iter != iterEnd; ++iter) { sd += pow( double(*iter - mean) , 2.0 ); }
  featureSd.push_back( sqrt(sd / double(length)) );

  // Calculate the median and quartiles. If the length is odd, this is the middle value.
  if (length % 2) {
    int middle = (length + 1) / 2 - 1;
    featureMedian.push_back( *(newCoverage.begin() + middle));

    // Now calculate the first and third quartiles. If there are (4n + 1) points,
    if (remainder(length - 1, 4) == 0) {
      int n = (length - 1) / 4;
      featureQ1.push_back((0.25 * *(newCoverage.begin() + n - 1)) + (0.75 * *(newCoverage.begin() + n)));
      featureQ3.push_back((0.75 * *(newCoverage.begin() + (3 * n))) + (0.25 * *(newCoverage.begin() + (3 * n) + 1)));
    } else if (remainder(length - 3, 4) == 0) {
      int n = (length - 3) / 4;
      featureQ1.push_back((0.75 * *(newCoverage.begin() + n)) + (0.25 * *(newCoverage.begin() + n + 1)));
      featureQ3.push_back((0.25 * *(newCoverage.begin() + (3 * n) + 1)) + (0.75 * *(newCoverage.begin() + (3 * n) + 2)));
    } else {
      cout << "Mathematical error in quartile calculation" << endl;
      exit(1);
    }

  // If the length is even, take the average of the middle two values.
  } else {
    int value1 = *(newCoverage.begin() + (length / 2) - 1);
    int value2 = *(newCoverage.begin() + ( (length / 2)) );
    featureMedian.push_back( (double(value1) + value2) / 2);

    // And determine the quartiles. If the lower half of the data has an odd length,
    // the quartiles are the midpoint of the lower and upper halves.
    if ( (length / 2) % 2) {
      int middle = ((length / 2) + 1) / 2;
      featureQ1.push_back( *(newCoverage.begin() + middle - 1));
      featureQ3.push_back( *(newCoverage.begin() + middle + (length / 2) ) - 1);

    // If the lower half has an even length, take the average of the middle two points.
    } else {
      int value1 = *(newCoverage.begin() + (length / 4) - 1);
      int value2 = *(newCoverage.begin() + (length / 4));
      featureQ1.push_back( (double(value1) + value2) / 2);
      value1 = *(newCoverage.begin() + 3 * (length / 4) - 1);
      value2 = *(newCoverage.begin() + 3 * (length / 4));
      featureQ3.push_back( (double(value1) + value2) / 2);
    }
  }

  // Now store the minimum and maximum values in the feature.
  featureMin.push_back(*(newCoverage.begin()));
  featureMax.push_back(*(newCoverage.end() - 1));
}

// Calculate the same values at the gene level.
void coverageData::processGene() {

  // Initialise variables.
  int length = geneCoverage.size();
  int sum    = 0;
  double sd  = 0.;

  // Sort the gene level coverage.
  std::sort(geneCoverage.begin(), geneCoverage.end());

  // Calculate values.
  vector<int>::iterator iter    = geneCoverage.begin();
  vector<int>::iterator iterEnd = geneCoverage.end();
  geneMin = *iter;
  geneMax = *(iterEnd - 1);

  // The mean for the gene.
  for (; iter != iterEnd; ++iter) { sum += *iter; }
  geneMean = sum / double(length);

  // The standard deviation for the gene.
  iter = geneCoverage.begin();
  for (; iter != iterEnd; ++iter) { sd += pow( double(*iter - geneMean) , 2.0 ); }
  geneSd = sqrt(sd / length);

  // The gene median.
  if (length % 2) {
    int middle = (length + 1) / 2 - 1;
    geneMedian = *(geneCoverage.begin() + middle);

    // Now calculate the first and third quartiles. If there are (4n + 1) points,
    if (remainder(length - 1, 4) == 0) {
      int n = (length - 1) / 4;
      geneQ1 = ((0.25 * *(geneCoverage.begin() + n - 1)) + (0.75 * *(geneCoverage.begin() + n)));
      geneQ3 = ((0.75 * *(geneCoverage.begin() + (3 * n))) + (0.25 * *(geneCoverage.begin() + (3 * n) + 1)));
    } else if (remainder(length - 3, 4) == 0) {
      int n = (length - 3) / 4;
      geneQ1 = ((0.75 * *(geneCoverage.begin() + n)) + (0.25 * *(geneCoverage.begin() + n + 1)));
      geneQ3 = ((0.25 * *(geneCoverage.begin() + (3 * n) + 1)) + (0.75 * *(geneCoverage.begin() + (3 * n) + 2)));
    } else {
      cout << "Mathematical error in quartile calculation" << endl;
      exit(1);
    }
  } else {
    int value1 = *(geneCoverage.begin() + (length / 2) - 1);
    int value2 = *(geneCoverage.begin() + ( (length / 2)) );
    geneMedian = (double(value1) + value2) / 2;

    // And determine the quartiles. If the lower half of the data has an odd length,
    // the quartiles are the midpoint of the lower and upper halves.
    if ( (length / 2) % 2) {
      int middle = ((length / 2) + 1) / 2;
      geneQ1 = ( *(geneCoverage.begin() + middle - 1));
      geneQ3 = ( *(geneCoverage.begin() + middle + (length / 2) ) - 1);

    // If the lower half has an even length, take the average of the middle two points.
    } else {
      int value1 = *(geneCoverage.begin() + (length / 4) - 1);
      int value2 = *(geneCoverage.begin() + (length / 4));
      geneQ1 = ( (double(value1) + value2) / 2);
      value1 = *(geneCoverage.begin() + 3 * (length / 4) - 1);
      value2 = *(geneCoverage.begin() + 3 * (length / 4));
      geneQ3 = ( (double(value1) + value2) / 2);
    }
  }
}
