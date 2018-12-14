//
//  fastqQualDiagnosing_lite.hpp
//  fastqQualDiagnosing_lite
//
//  Created by Hao Yu on 6/13/16.
//  Copyright © 2016 Hao Yu. All rights reserved.
//

#ifndef fastqQualDiagnosing_lite_hpp
#define fastqQualDiagnosing_lite_hpp

#include <string>
#include <map>

#endif /* fastqQualDiagnosing_lite_hpp */

using namespace std;

typedef struct {
    unsigned long int all_reads;
    unsigned long int all_bases;
         unsigned int read_len;
    
    unsigned long int qual_33_bases;
    unsigned long int qual_33_64_bases;
    unsigned long int qual_64_bases;
    unsigned long int qual_tooHigh_bases;
    
    float qual_tooLow_ratio;
    float qual_33_ratio;
    float qual_33_64_ratio;
    float qual_64_ratio;
    float qual_tooHigh_ratio;
} fastq_box;

fastq_box statistics_fastq(string line);
