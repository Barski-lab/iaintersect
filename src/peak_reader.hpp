#ifndef PEAK_READER_HPP
#define PEAK_READER_HPP

#include "config.hpp"
#include "string_tools.hpp"

class GeneInfo {
    std::string refseq_id;
    std::string gene_id;
    long txStart;
    long txEnd;
    bool strand;
    std::string region;
    GeneInfo (std::string new_refseq_id, std::string new_gene_id, long new_txStart, long new_txEnd, bool new_strand, std::string new_region)
        : refseq_id (new_refseq_id)
        , gene_id (new_gene_id)
        , txStart (new_txStart)
        , txEnd (new_txEnd)
        , strand (new_strand)
        , region (new_region)
    {

    }
};

class PeakRecord {
    std::string chr;
    long start;
    long end;
    int length;
    long abs_summit;
    double pileup;
    double log10p;
    double fold_enrichment;
    double log10q;
    std::string name;
    GeneInfo gene_info;
    PeakRecord (const std::string & line);
};

typedef boost::shared_ptr<PeakRecord> PeakRecordPtr;

struct Coord{
    std::string chr;
    long start;
    long end;
    Coord (std::string new_chr, long new_start, long new_end)
        : chr (new_chr)
        , start (new_start)
        , end (new_end)
    {
    }

};

class PeakReader {
private:
    // [chromosome][peak start][peak end] = vector of PeakRecord
    std::map<std::string, std::map<long, std::map<long, std::vector<PeakRecordPtr> > > > peak_data;
    bool include_chr(const std::string chr);
    bool include_start(const std::string chr, long start);
    bool include_end(const std::string chr, long start, long end);
    std::map<std::string, std::map<long, std::map<long, std::vector<PeakRecordPtr> > > >::iterator chr_it;
    std::map<long, std::map<long, std::vector<PeakRecordPtr> > > start_it;
    std::map<long, std::vector<PeakRecordPtr> > end_it;
public:
    bool load (const std::string & full_filename);
    bool next();
    std::vector<PeakRecordPtr> value();
    Coord coordinates();
    void reset(); // reset iterators and current_record for next function
    bool update (std::string chr, long start, long end, GeneInfo new_gene_info);
    PeakReader (const std::string & full_filename);
}


#endif // PEAK_READER_HPP
