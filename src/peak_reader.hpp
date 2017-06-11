#ifndef PEAK_READER_HPP
#define PEAK_READER_HPP

#include "config.hpp"
#include "string_tools.hpp"

class GeneInfo {
public:
    std::string refseq_id;
    std::string gene_id;
    int txStart;
    int txEnd;
    std::string strand;
    std::string region;

    GeneInfo (std::string new_refseq_id, std::string new_gene_id, int new_txStart, int new_txEnd, std::string new_strand, std::string new_region)
        : refseq_id (new_refseq_id)
        , gene_id (new_gene_id)
        , txStart (new_txStart)
        , txEnd (new_txEnd)
        , strand (new_strand)
        , region (new_region)
    {
    };

    GeneInfo()
        : refseq_id ("NULL")
        , gene_id ("NULL")
        , txStart (0)
        , txEnd (0)
        , strand ("NULL")
        , region ("NULL")
    {

    };
    void print ();
};

class PeakRecord {
public:
    std::string chr;
    int start;
    int end;
    int length;
    int abs_summit;
    double pileup;
    double log10p;
    double fold_enrichment;
    double log10q;
    std::string name;
    GeneInfo gene_info;
    bool broad_flag;
    PeakRecord (const std::string & line, bool broad);
    void print();
};

typedef boost::shared_ptr<PeakRecord> PeakRecordPtr;

struct Coord{
    std::string chr;
    int start;
    int end;
    Coord (std::string new_chr, int new_start, int new_end)
        : chr (new_chr)
        , start (new_start)
        , end (new_end)
    {
    }

};

class PeakReader {
private:
    // [chromosome][peak start][peak end] = vector of PeakRecord
    std::map<std::string, std::map<int, std::map<int, std::vector<PeakRecordPtr> > > > peak_data;
    bool include_chr(const std::string chr);
    bool include_start(const std::string chr, int start);
    bool include_end(const std::string chr, int start, int end);
    std::map<std::string, std::map<int, std::map<int, std::vector<PeakRecordPtr> > > >::iterator chr_it;
    std::map<int, std::map<int, std::vector<PeakRecordPtr> > >::iterator start_it;
    std::map<int, std::vector<PeakRecordPtr> >::iterator end_it;
    std::string filename;
    bool broad_flag;
public:
    bool load (const std::string & full_filename);
    bool save (const std::string & output_filename);
    bool next();
    bool has_next();
    std::vector<PeakRecordPtr> value();
    Coord coordinates();
    void print (ostream& output_stream);
    std::string get_filename();
    void reset(); // reset iterators and current_record for next function
    bool update (std::string chr, int start, int end, const GeneInfo & new_gene_info);
    PeakReader (const std::string & full_filename);
};


#endif // PEAK_READER_HPP
