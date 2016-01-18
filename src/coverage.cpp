#include "api/BamMultiReader.h"
#include "dataProcessing.h"
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
using namespace BamTools;

// Parse a file and add all lines to the list of regions.
vector<string> getRegions(string file) {
  vector<string> list;
  string line;
  ifstream infile(file.c_str());;

  // Open the file reading.
  while (getline(infile, line)) {list.push_back(line);}

  return list;
}

// Check that the region string is valid.
bool ParseRegionString(const string& regionString, const BamMultiReader& reader, BamRegion& region) {

  // Check first for empty string.
  if ( regionString.empty() ) {
    return false;
  }
  
  // Look for a colon.
  size_t foundFirstColon = regionString.find(':');
  
  // Store chrom strings, and numeric positions.
  string startChrom;
  string stopChrom;
  int startPos;
  int stopPos;
  
  // If no colon is found, use entire contents of requested chromosome.
  // Store the entire region string as startChrom name and use BamReader methods
  // to check if its valid for current BAM file.
  if ( foundFirstColon == string::npos ) {
    startChrom = regionString;
    startPos   = 0;
    stopChrom  = regionString;
    stopPos    = -1;
  }
  
  // If a colon is found, we have some sort of startPos requested.
  else {
    
    // Store the start chrom from beginning to first colon.
    startChrom = regionString.substr(0, foundFirstColon);
    
    // Look for ".." or "-"  after the colon.
    size_t foundRangeDots = regionString.find("..", foundFirstColon + 1);
    size_t foundDash      = regionString.find("-", foundFirstColon + 1);
    
    // If no dots or dash found, we have a startPos but no range. Store the 
    // contents before colon as startChrom, after as startPos.
    if ( foundRangeDots == string::npos and foundDash == string::npos ) {
      startPos   = atoi( regionString.substr(foundFirstColon + 1).c_str() ); 
      stopChrom  = startChrom;
      stopPos    = -1;
    } 
    
    // A ".." or "-" is found, so we have some sort of range selected.
    else {
      size_t foundSecondColon;
      size_t foundRange;
      size_t offset;
      if ( foundRangeDots == string::npos) {
        foundRange = foundDash;
        offset     = 1;
      }
      else { 
        foundRange = foundRangeDots;
        offset     = 2;
      }

      // Store the startPos between first colon and range dots ".." or "-" and look for a
      // second colon.
      startPos         = atoi( regionString.substr(foundFirstColon + 1, foundRange - foundFirstColon - 1).c_str() );
      foundSecondColon = regionString.find(':', foundRange + 1);
      
      // If no second colon found, so we have a "standard" chrom:start..stop input format (on single chrom).
      if ( foundSecondColon == string::npos ) {
        stopChrom  = startChrom;
        stopPos    = atoi( regionString.substr(foundRange + offset).c_str() );
      }
      
      // If a second colon is found, we have a range requested across 2 chrom's.
      else {
        stopChrom  = regionString.substr(foundRange + offset, foundSecondColon - ( foundRange + offset ));
        stopPos    = atoi( regionString.substr(foundSecondColon + 1).c_str() );
      }
    }
  }

  // Validate reference IDs & genomic positions.
  const RefVector references = reader.GetReferenceData();

  // If startRefID not found, return false.
  int startRefID = reader.GetReferenceID(startChrom);
  if ( startRefID == -1 ) return false;

  // startPos cannot be greater than or equal to reference length.
  const RefData& startReference = references.at(startRefID);
  if ( startPos >= startReference.RefLength ) return false;

  // If stopRefID not found, return false.
  int stopRefID = reader.GetReferenceID(stopChrom);
  if ( stopRefID == -1 ) return false;

  // stopPosition cannot be larger than reference length.
  const RefData& stopReference = references.at(stopRefID);
  if ( stopPos > stopReference.RefLength ) return false;

  // If no stopPosition specified, set to reference end.
  if ( stopPos == -1 ) stopPos = stopReference.RefLength;

  // Set up the Region struct & return.
  region.LeftRefID     = startRefID;
  region.LeftPosition  = startPos;
  region.RightRefID    = stopRefID;;
  region.RightPosition = stopPos;
  return true;
}

bool processCigar(BamAlignment& al, BamRegion& region, int startPosition, vector<int>& coverage) {
  
  // Intialize local variables.
  const int numCigarOps = (const int)al.CigarData.size(); 
  int positionInRegion  = al.Position - startPosition;
  
  // Iterate over the CIGAR operations. 
  for (int i = 0; i < numCigarOps; ++i ) {
    const CigarOp& op = al.CigarData.at(i);
  
    // If the CIGAR string indicates a match.
    if (op.Type == 'M') {
      for (int j = positionInRegion; j < (positionInRegion + (int)op.Length); ++j) {
        if (j < coverage.size()) {coverage[j]++;}
      }
      positionInRegion += (int)op.Length;
    }
  
    // If the bases are soft or hard clipped bases, advance the position in the region, but do not update
    // coverage.
    else if (op.Type == 'S' || op.Type == 'H') { positionInRegion += (int)op.Length; }
  
    // If there is an insertion, do nothing. All the bases in the insertion do not cover reference bases, so 
    // should not be counted and the position in the region is not advanced.
    else if (op.Type == 'I') { }
  
    // If there is an deletion, count the deleted bases as covered and advance the position in the region.
    else if (op.Type == 'D') {
      for (int j = positionInRegion; j < (positionInRegion + (int)op.Length); ++j) {
        if (j < coverage.size()) {coverage[j]++;}
      }
      positionInRegion += (int)op.Length;
    }
  }
  return true;
}

int main(int argc, char * argv[])
{
  // record command line parameters
  string commandLine = argv[0];

  // Required variables.
  int c;
  string gene = "";
  string transcript;
  string regionsFile;
  string output;
  vector<string> inputFiles;

  static struct option long_options[] =
    {
      {"help", no_argument, 0, 'h'},
      {"bam", required_argument, 0, 'b'},
      {"regions", required_argument, 0, 'r'},
      {"output", required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };

  while (true) {
    int option_index = 0;
    c = getopt_long(argc, argv, "hb:g:t:r:o:", long_options, &option_index);

    if (c == -1) // end of options
      break;

    switch (c) {

      // The BAM file..
      case 'b':
        inputFiles.push_back(optarg);
        break;

      // The list of regions.
      case 'r':
        regionsFile = optarg;
        break;

      // The output file.
      case 'o':
        output = optarg;
        break;

      default:
        abort ();
    }
  }

  // Bam files must be specified.
  if (inputFiles.size() == 0) {
    cerr << "Please specify a BAM file or files (--bam, -b)." << endl;
    exit(1);
  }

  // A file containing a list of regions must be specified.
  if (regionsFile == "" ) {
    cerr << "Please specify a file containing a list of regions (--regions, -r)." << endl;
    exit(1);
  }

  // Read the file containing regions and add all regions to the list.
  vector<string> regionList = getRegions(regionsFile);

  // Define a structure for holding mean information.
  coverageData cov(regionList.size());

  // Open the multireader
  BamMultiReader reader;
  if ( !reader.Open(inputFiles) ) {
    cerr << "bamtools count ERROR: could not open input BAM file(s)... Aborting." << endl;
    exit(1);
  }

  // Define a region and alignment.
  BamRegion region;

  // Retrieve references.
  BamTools::RefVector references = reader.GetReferenceData();

  // Open output file (or stdout) for writing.
  ofstream outputFile;
  bool outputRequested = false;
  if (output != "") {
    outputRequested = true;
    outputFile.open(output.c_str(), std::ios::out);
  }
  ostream & outFile = ( outputRequested ? outputFile : cout);

  // Create a vector to store the info for each transcript.
  vector<double> means;
  vector<double> sds;

  // Variables for defining a feature id.
  int exonId = 1;

  // Loop over all regions.
  vector<string>::iterator iter    = regionList.begin();
  vector<string>::iterator iterEnd = regionList.end();
  for (; iter != iterEnd; ++iter) {

    // Check that the region is valid and locate indexes.
    if ( !ParseRegionString(*iter, reader, region) ) {
      cerr << "ERROR: Invalid region string: " << *iter << endl;
      exit(1);
    }

    // Attempt to find index files.
    reader.LocateIndexes();

    // If index data available for all BAM files, we can use SetRegion.
    if ( reader.HasIndexes() ) {

      // Attempt to set region on reader.
      if ( !reader.SetRegion(region.LeftRefID, region.LeftPosition, region.RightRefID, region.RightPosition) ) {
        cerr << "bamtools count ERROR: set region failed. Check that REGION describes a valid range" << endl;
        reader.Close();
        exit(1);
      }

      // Determine the length of the region and use this to define the start of each feature in the array of
      // coverage data.
      unsigned int length = region.RightPosition - region.LeftPosition + 1;
      cov.featureLengths.push_back(length);

      // Create an id with the region included.
      ostringstream oss;
      oss << exonId << "\t" << *iter;
      cov.ids.push_back(oss.str());
      exonId++;

      // Set up the pileup engine.
      //PileupEngine pileup;

      // Parse the region and build up statistics.
      cout << "REGION: " << region.LeftPosition << "-" << region.RightPosition << endl;

      // Define a new BamAlignment. Declaring here will ensure that if this region has no reads, but the previous
      // region did, the alignment object will be cleared.
      BamAlignment al;

      // Get the first alignment to set the start coordinate of the first base in the first read.
      reader.GetNextAlignment(al);

      // If there are no alignments.
      if (al.Position == -1) {
        cout << "  NO ALIGNMENTS -> NO COVERAGE" << endl;
      } else {

        // Initialise variables.
        vector<int> coverage(region.RightPosition - al.Position);
        int coverageStart = al.Position;

        // Process the first read.
        processCigar(al, region, coverageStart, coverage);

        // Loop over the remaining reads spanning the region.
        while ( reader.GetNextAlignment(al) ) { processCigar(al, region, coverageStart, coverage); }

        //
        int start = region.LeftPosition - coverageStart;
        vector<int>::iterator covIter = coverage.begin() + start;
        vector<int>::iterator covIterEnd = coverage.end();
      }
    }
  }

  // Sort the coverages in each feature (for median calculation) and generate a sorted list of all depths (for
  // the gene level calculation).
  //cov.sort();

  // Iterators for features.
  //vector<string>::iterator idIter    = cov.ids.begin();
  //vector<string>::iterator idIterEnd = cov.ids.end();

  //vector<int>::iterator minIter    = cov.featureMin.begin();
  //vector<int>::iterator minIterEnd = cov.featureMin.end();

  //vector<int>::iterator maxIter    = cov.featureMax.begin();
  //vector<int>::iterator maxIterEnd = cov.featureMax.end();

  //vector<double>::iterator meanIter    = cov.featureMean.begin();
  //vector<double>::iterator meanIterEnd = cov.featureMean.end();

  //vector<double>::iterator medIter    = cov.featureMedian.begin();
  //vector<double>::iterator medIterEnd = cov.featureMedian.end();

  //vector<double>::iterator sdIter    = cov.featureSd.begin();
  //vector<double>::iterator sdIterEnd = cov.featureSd.end();

  // Include a header line.
  //outFile << "#id\tregion\tmin\tmax\tmedian\tmean\tsd" << endl;

  // Iterate over the feature minimum values and increment all other iterators as we go.
  //for (; idIter != idIterEnd; ++idIter) {
    //outFile << *idIter << "\t" << *minIter << "\t" << *maxIter << "\t" << *medIter << "\t" << *meanIter << "\t" << *sdIter << endl;

    // Increment the exon id.
    //exonId++;

    // Increment the iterators.
    //++minIter;
    //++maxIter;
    //++meanIter;
    //++medIter;
    //++sdIter;
  //}

  // Now include the gene level information.
  //outFile << "gene\tNA\t" << cov.geneMin << "\t" << cov.geneMax << "\t" << cov.geneMedian << "\t" << cov.geneMean << "\t" << cov.geneSd << endl;
}
