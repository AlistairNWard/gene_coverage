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
  depthData.clear();
}

coverageData::~coverageData(void) {
}

void coverageData::update(int value) {

  // Add the coverage value to the end of the array.
  depthData.push_back(value);
}

// Sort the data according to feature.
void coverageData::sort(void) {

  // Define iterators.
  vector<int>::iterator iiter    = featureLengths.begin();
  vector<int>::iterator iiterEnd = featureLengths.end();

  vector<int>::iterator depthIter;
  vector<int>::iterator depthIterEnd;

  int length    = 0;
  int start     = 0;
  int end       = 0;
  int sum       = 0;
  double sdTemp = 0.;
  cout << "DATA: " << featureLengths.size() << " " << depthData.size() << endl;
  for (; iiter != iiterEnd; ++iiter) {

    // Define the end of the feature and the length.
    end    = start + *iiter;
    length = end - start;
    std::sort(depthData.begin() + start, depthData.begin() + end);

    // Define iterators for feature.
    depthIter    = depthData.begin() + start;
    depthIterEnd = depthData.begin() + end;

    // Find the min and max values
    featureMin.push_back(*depthIter);
    featureMax.push_back(*(depthIterEnd - 1));

    cout << "    " << end << " " << length << " " << *depthIter << " " << *(depthIterEnd - 1) << " " << depthData.size() << endl;

    // Calculate the mean.
    sum = 0;
    for (; depthIter != depthIterEnd; ++depthIter) { sum += *depthIter; }
    double mean = sum / double(length);
    featureMean.push_back(mean);

    // and the standard deviation.
    sdTemp    = 0;
    depthIter = depthData.begin() + start;
    for (; depthIter != depthIterEnd; ++depthIter) { sdTemp += pow( double(*depthIter - mean) , 2.0 ); }
    featureSd.push_back( sqrt(sdTemp / double(length)) );

    // Calculate the median. If the length is odd, this is the middle value.
    if (length % 2) { featureMedian.push_back( *(depthData.begin() + start + ((length + 1) / 2)) );

    // If the length is even, take the average of the middle two values.
    } else {
      int value1 = *(depthData.begin() + start + (length / 2));
      int value2 = *(depthData.begin() + start + ( (length / 2) - 1) );
      featureMedian.push_back( (double(value1) + value2) / 2);
    }

    // Update the start.
    start = end;
  }

  // Calculate the same values for the gene.
  geneMin = *min_element(featureMin.begin(), featureMin.end());
  geneMax = *max_element(featureMax.begin(), featureMax.end());

  // The mean for the gene.
  vector<double>::iterator meanIter    = featureMean.begin();
  vector<double>::iterator meanIterEnd = featureMean.end();
  double meanSum = 0.;
  for (; meanIter != meanIterEnd; ++meanIter) { meanSum += *meanIter; }
  geneMean = meanSum / featureMean.size();

  // The standard deviation for the gene.
  depthIter    = depthData.begin();
  depthIterEnd = depthData.end();
  sdTemp       = 0.;
  for (; depthIter != depthIterEnd; ++depthIter) { sdTemp += pow( double(*depthIter - geneMean) , 2.0 ); }
  geneSd = sqrt(sdTemp / depthData.size());

  // The gene median.
  std::sort(depthData.begin(), depthData.end());
  depthIter = depthData.begin();
  length    = depthData.size();
  if (length % 2) { geneMedian = *(depthData.begin() + ( (length + 1) / 2 ));
  } else {
    int value1 = *(depthData.begin() + (length / 2));
    int value2 = *(depthData.begin() + ( (length / 2) -1) );
    geneMedian = (double(value1) + value2) / 2;
  }
}
