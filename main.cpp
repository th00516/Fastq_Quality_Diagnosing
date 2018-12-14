//
//  main.cpp
//  fastqQualDiagnosing_lite
//
//  Created by Hao Yu on 6/13/16.
//  Copyright © 2016 Hao Yu. All rights reserved.
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "fastqQualDiagnosing_lite.hpp"

int main(int argc, const char *argv[])
{
    int sample_line;
    if (argv[2])
    {
        sample_line = atoi(argv[2]) * 4;
    }
    else
    {
        sample_line = 100000 * 4;
    }
    
    fastq_box fastq;
    
    string path;
    string filename;

    fstream fs;
    fs.open(argv[1], fstream::in);
    if (fs.fail())
    {
        cerr << "==== fastqQualDiagnosing_lite ===="                                   << '\n'
             << "Version 1.1, Created by:"                                             << '\n'
             << setw(8) << "" << "Hao Yu  (yuhao@genomics.cn)"                         << '\n'
             << setw(8) << "" << "Qiye Li (liqiye@genomics.cn)"                        << '\n'
                                                                                       << '\n'
             << setw(2) << "USAGE: " << argv[0] << " <fastq_file> [sample_head_lines]" << '\n'
             << endl;
        
        return 1;
    }
    else
    {
        path = argv[1];
        filename = path.substr(path.find_last_of('/') + 1);
    }
    
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(fs);
    
    unsigned long int nr = 1;
    string line;
    while (getline(in, line))
    {
        if (nr > sample_line)
        {
            break;
        }
        
        if (nr % 4 == 0)
        {
            // statistics the information of fastq file
            fastq = statistics_fastq(line);
        }
        
        nr ++;
    }
    
    fs.close();
    
    fastq.qual_33_ratio      = float(fastq.qual_33_bases)      / float(fastq.all_bases);
    fastq.qual_33_64_ratio   = float(fastq.qual_33_64_bases)   / float(fastq.all_bases);
    fastq.qual_64_ratio      = float(fastq.qual_64_bases)      / float(fastq.all_bases);
    fastq.qual_tooHigh_ratio = float(fastq.qual_tooHigh_bases) / float(fastq.all_bases);
    
    fastq.all_reads = nr / 4;
    
    // output result
    cout << "================================================================" << '\n';
    cout << "FASTQ file name: " << filename                                    << '\n';
    cout << "================================================================" << '\n';
    cout << endl;
    
    cout << right << setw(40) << "Summary: Total number of sampled reads: " << setw(12) << fastq.all_reads << '\n'
                  << setw(40) << "Total number of sampled bases: "          << setw(12) << fastq.all_bases << '\n'
         << endl;
    if (fastq.all_bases % fastq.all_reads == 0)
    {
        fastq.read_len = int(fastq.all_bases) / int(fastq.all_reads);
        cout << right << setw(40) << "Length of each read: " << setw(12) << fastq.read_len  << '\n';
        cout << endl;
    }
    else
    {
        cout << "Warning: Not all the reads are the same length." << '\n'
             << setw(9) << ""
             << "Your FASTQ file may be a mixture of FASTQ files from separate source." << '\n';
        cout << endl;
    }
    
    
    cout << right << setw(10) << "ASCII_RANGE"
                  << setw(15) << "BASE_NUM"
                  << setw(10) << "RATIO"
                  << setw(4)  << ""
         << left  << setw(35) << "REMARKS"
         << '\n';
    
    cout << right << setw(12) << "[33 , 63]:"
                  << setw(14)                    << fastq.qual_33_bases
         << fixed << setw(10) << setprecision(4) << fastq.qual_33_ratio
                  << setw(4)  << ""
         << left  << setw(35) << "Illumina 1.8+ Phred+33"
         << '\n';
    
    cout << right << setw(12) << "[64 , 74]:"
                  << setw(14)                    << fastq.qual_33_64_bases
         << fixed << setw(10) << setprecision(4) << fastq.qual_33_64_ratio
                  << setw(4)  << ""
         << left  << setw(35) << "Illumina 1.8+ Phred+33/1.3+ Phred+64"
         << '\n';
    
    cout << right << setw(12) << "[75 ,105]:"
                  << setw(14)                    << fastq.qual_64_bases
         << fixed << setw(10) << setprecision(4) << fastq.qual_64_ratio
                  << setw(4)  << ""
         << left  << setw(35) << "Illumina 1.3+ Phred+64"
         << '\n';
    
    cout << right << setw(12) << "[106,...):"
                  << setw(14)                    << fastq.qual_tooHigh_bases
         << fixed << setw(10) << setprecision(4) << fastq.qual_tooHigh_ratio
         << '\n';
    
    cout << endl;
    
    // output the conclusion
    if      (fastq.qual_33_ratio > fastq.qual_64_ratio)
    {
        if      (fastq.qual_33_ratio >= 0.3 && fastq.qual_64_ratio <  0.0001)
        {
            cout << "======== Conclusion ========"  << '\n';
            cout << "These reads may be set quality based on Illumina 1.8+ Phred+33." << '\n';
            cout << endl;
        }
        else if (fastq.qual_33_ratio >= 0.3 && fastq.qual_64_ratio >= 0.0001)
        {
            cout << "======== Warning ========" << '\n';
            cout << "Warning: Your Illumina 1.8+ Phred+33 reads may be mixed up with Illumina 1.3+ Phred+64 reads." << '\n';
            cout << endl;
        }
        else if (fastq.qual_33_ratio <  0.3 && fastq.qual_64_ratio >= 0.0001)
        {
            cout << "======== Fatal Error ========" << '\n';
            cout << "Warning: Your Illumina 1.8+ Phred+33 reads may be mixed up with too many Illumina 1.3+ Phred+64 reads."  << '\n'
                 << setw(9) << ""
                 << "I suggest you to cautiously check your FASTQ file, if too hard, you should abandon this doubtable data." << '\n';
            cout << endl;
        }
    }
    else if (fastq.qual_64_ratio > fastq.qual_33_ratio)
    {
        if      (fastq.qual_64_ratio >= 0.3 && fastq.qual_33_ratio <  0.0001)
        {
            cout << "======== Conclusion ========"  << '\n';
            cout << "These reads may be set quality based on Illumina 1.3+ Phred+64." << '\n';
            cout << endl;
        }
        else if (fastq.qual_64_ratio >= 0.3 && fastq.qual_33_ratio >= 0.0001)
        {
            cout << "======== Warning ========" << '\n';
            cout << "Warning: Your Illumina 1.3+ Phred+64 reads may be mixed up with Illumina 1.8+ Phred+33 reads." << '\n';
            cout << endl;
        }
        else if (fastq.qual_64_ratio <  0.3 && fastq.qual_33_ratio >= 0.0001)
        {
            cout << "======== Fatal Error ========" << '\n';
            cout << "Warning: Your Illumina 1.3+ Phred+64 reads may be mixed up with too many Illumina 1.8+ Phred+33 reads."  << '\n'
                 << setw(9) << ""
                 << "I suggest you to cautiously check your FASTQ file, if too hard, you should abandon this doubtable data." << '\n';
            cout << endl;
        }
    }
    
    if (fastq.qual_33_64_ratio > 0.5)
    {
        cout << "======== Warning ========" << '\n';
        cout << "Too many bases are too ambiguous in your FASTQ file (QUAL_ASCII: [64 , 73])," << '\n'
             << "this ambiguousness may obstruct your judgment." << '\n';
        cout << endl;
    }
    
    return 0;
}
