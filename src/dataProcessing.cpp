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
  //depthData.clear();
}

coverageData::~coverageData(void) {
}

// If a feature has no coverage, add the statistics to the correct fields.
void coverageData::noCoverage() {
  featureMin.push_back(0);
  featureMax.push_back(0);
  featureMean.push_back(0);
  featureMedian.push_back(0);
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
  vector<int>::iterator iter    = coverage.begin() + start;
  vector<int>::iterator iterEnd = coverage.end();
  for (; iter != iterEnd; ++iter) {

    // Keep a sum of the coverage values.
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

  // Calculate the median. If the length is odd, this is the middle value.
  if (length % 2) { featureMedian.push_back( *(newCoverage.begin() + ((length + 1) / 2)) );

  // If the length is even, take the average of the middle two values.
  } else {
    int value1 = *(newCoverage.begin() + (length / 2));
    int value2 = *(newCoverage.begin() + ( (length / 2) - 1) );
    featureMedian.push_back( (double(value1) + value2) / 2);
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
  if (length % 2) { geneMedian = *(geneCoverage.begin() + ( (length + 1) / 2 ));
  } else {
    int value1 = *(geneCoverage.begin() + (length / 2));
    int value2 = *(geneCoverage.begin() + ( (length / 2) -1) );
    geneMedian = (double(value1) + value2) / 2;
  }
}




//void coverageData::update(int value) {
//
//  // Add the coverage value to the end of the array.
//  depthData.push_back(value);
//}

// Sort the data according to feature.
//void coverageData::sort(void) {
//
//  // Define iterators.
//  vector<int>::iterator iiter    = featureLengths.begin();
//  vector<int>::iterator iiterEnd = featureLengths.end();
//
//  vector<int>::iterator depthIter;
//  vector<int>::iterator depthIterEnd;
//
//  int length    = 0;
//  int start     = 0;
//  int end       = 0;
//  int sum       = 0;
//  double sdTemp = 0.;
//  cout << "DATA: " << featureLengths.size() << " " << depthData.size() << endl;
//  for (; iiter != iiterEnd; ++iiter) {
//
//    // Define the end of the feature and the length.
//    end    = start + *iiter;
//    length = end - start;
//    std::sort(depthData.begin() + start, depthData.begin() + end);
//
//    // Define iterators for feature.
//    depthIter    = depthData.begin() + start;
//    depthIterEnd = depthData.begin() + end;
//
//    // Find the min and max values
//    featureMin.push_back(*depthIter);
//    featureMax.push_back(*(depthIterEnd - 1));
//
//    cout << "    " << end << " " << length << " " << *depthIter << " " << *(depthIterEnd - 1) << " " << depthData.size() << endl;
//
//    // Calculate the mean.
//    sum = 0;
//    for (; depthIter != depthIterEnd; ++depthIter) { sum += *depthIter; }
//    double mean = sum / double(length);
//    featureMean.push_back(mean);
//
//    // and the standard deviation.
//    sdTemp    = 0;
//    depthIter = depthData.begin() + start;
//    for (; depthIter != depthIterEnd; ++depthIter) { sdTemp += pow( double(*depthIter - mean) , 2.0 ); }
//    featureSd.push_back( sqrt(sdTemp / double(length)) );
//
//    // Calculate the median. If the length is odd, this is the middle value.
//    if (length % 2) { featureMedian.push_back( *(depthData.begin() + start + ((length + 1) / 2)) );
//
//    // If the length is even, take the average of the middle two values.
//    } else {
//      int value1 = *(depthData.begin() + start + (length / 2));
//      int value2 = *(depthData.begin() + start + ( (length / 2) - 1) );
//      featureMedian.push_back( (double(value1) + value2) / 2);
//    }
//
//    // Update the start.
//    start = end;
//  }
//
//  // Calculate the same values for the gene.
//  geneMin = *min_element(featureMin.begin(), featureMin.end());
//  geneMax = *max_element(featureMax.begin(), featureMax.end());
//
//  // The mean for the gene.
//  vector<double>::iterator meanIter    = featureMean.begin();
//  vector<double>::iterator meanIterEnd = featureMean.end();
//  double meanSum = 0.;
//  for (; meanIter != meanIterEnd; ++meanIter) { meanSum += *meanIter; }
//  geneMean = meanSum / featureMean.size();
//
//  // The standard deviation for the gene.
//  depthIter    = depthData.begin();
//  depthIterEnd = depthData.end();
//  sdTemp       = 0.;
//  for (; depthIter != depthIterEnd; ++depthIter) { sdTemp += pow( double(*depthIter - geneMean) , 2.0 ); }
//  geneSd = sqrt(sdTemp / depthData.size());
//
//  // The gene median.
//  std::sort(depthData.begin(), depthData.end());
//  depthIter = depthData.begin();
//  length    = depthData.size();
//  if (length % 2) { geneMedian = *(depthData.begin() + ( (length + 1) / 2 ));
//  } else {
//    int value1 = *(depthData.begin() + (length / 2));
//    int value2 = *(depthData.begin() + ( (length / 2) -1) );
//    geneMedian = (double(value1) + value2) / 2;
//  }
//}
