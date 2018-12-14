//
//  fastqQualDiagnosing_lite.cpp
//  fastqQualDiagnosing_lite
//
//  Created by Hao Yu on 6/13/16.
//  Copyright © 2016 Hao Yu. All rights reserved.
//

#include <iostream>
#include "fastqQualDiagnosing_lite.hpp"

fastq_box fastq;

fastq_box statistics_fastq(string line)
{
    for (int site = 0; site < line.size(); site ++)
    {
        if      (line[site] >= '!' && line[site] <= '?')
        {
            fastq.qual_33_bases ++;
        }
        else if (line[site] >= '@' && line[site] <= 'J')
        {
            fastq.qual_33_64_bases ++;
        }
        else if (line[site] >= 'K' && line[site] <= 'i')
        {
            fastq.qual_64_bases ++;
        }
        else if (line[site] >= 'j')
        {
            fastq.qual_tooHigh_bases ++;
        }
    }
    fastq.all_bases += line.size();
        
    return fastq;
}
