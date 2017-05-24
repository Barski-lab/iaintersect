#include "peak_reader.hpp"

using namespace string_tools;

PeakRecord::PeakRecord(const std::string & line){

// chr2L<->18206<->18433<->228<--->18302<->52.00<->17.92226<------>4.15821>16.65316<------>SRR1198790.sorted_peak_1

    vector<string> line_splitted = split_line(line);
    if (line_splitted.size() < 9){
        throw ("PeakRecord constructor line length error");
    }
    chr = line_splitted[0];
    // start
    if (not str_to_int(start, line_splitted[1])){
        throw ("PeakRecord constructor error");
    }
    // end
    if (not str_to_int(end, line_splitted[2])){
        throw ("PeakRecord constructor error");
    }
    // length
    if (not str_to_int(length, line_splitted[3])){
        throw ("PeakRecord constructor error");
    }
    // abs_summit
    if (not str_to_int(abs_summit, line_splitted[4])){
        throw ("PeakRecord constructor error");
    }
    // pileup
    if (not str_to_double(pileup, line_splitted[5])){
        throw ("PeakRecord constructor error");
    }
    // log10p
    if (not str_to_double(log10p, line_splitted[6])){
        throw ("PeakRecord constructor error");
    }
    // fold_enrichment
    if (not str_to_double(fold_enrichment, line_splitted[7])){
        throw ("PeakRecord constructor error");
    }
    // log10q
    if (not str_to_double(log10q, line_splitted[8])){
        throw ("PeakRecord constructor error");
    }
    name = line_splitted[9];
}


PeakReader::PeakReader(const std::string & full_filename){
    if (!load (full_filename)){
        throw "PeakReader construcctor error: couldn't load any data";
    }
}


bool PeakReader::include_chr(const std::string chr){
    std::map<std::string, std::map<int, std::map<int, std::vector<PeakRecordPtr> > > >::iterator it = peak_data.find(chr);
    if(it != peak_data.end()) {
        return true;
    }
    return false;
}

bool PeakReader::include_start(const std::string chr, int start){
    if (include_chr(chr)){
        std::map<int, std::map<int, std::vector<PeakRecordPtr> > >::iterator it = peak_data[chr].find(start);
        if(it != peak_data[chr].end()) {
            return true;
        }
    }
    return false;
}

bool PeakReader::include_end(const std::string chr, int start, int end){
    if ( include_start(chr,start) ){
        std::map<int, std::vector<PeakRecordPtr> >::iterator it = peak_data[chr][start].find(end);
        if(it != peak_data[chr][start].end()) {
            return true;
        }
    }
    return false;
}

void PeakReader::reset(){
    chr_it = peak_data.begin();
    start_it = chr_it->second.begin();
    end_it = start_it->second.begin();
}

bool PeakReader::next(){
    std::map<std::string, std::map<int, std::map<int, std::vector<PeakRecordPtr> > > >::iterator chr_it_copy = chr_it;
    std::map<int, std::map<int, std::vector<PeakRecordPtr> > >::iterator start_it_copy = start_it;
    std::map<int, std::vector<PeakRecordPtr> >::iterator end_it_copy = end_it;
    if (++end_it_copy != start_it->second.end()){
        end_it++;
    } else if (++start_it_copy != chr_it->second.end()){
        start_it++;
        end_it = start_it->second.begin();
    } else if (++chr_it_copy != peak_data.end()){
        chr_it++;
        start_it = chr_it->second.begin();
        end_it = start_it->second.begin();
    } else {
        return false;
    }
    return true;
}

bool PeakReader::has_next(){
    std::map<std::string, std::map<int, std::map<int, std::vector<PeakRecordPtr> > > >::iterator chr_it_copy = chr_it;
    std::map<int, std::map<int, std::vector<PeakRecordPtr> > >::iterator start_it_copy = start_it;
    std::map<int, std::vector<PeakRecordPtr> >::iterator end_it_copy = end_it;
    if (++end_it_copy != start_it->second.end() or ++start_it_copy != chr_it->second.end() or ++chr_it_copy != peak_data.end()){
        return true;
    };
    return false;
}


std::vector<PeakRecordPtr> PeakReader::value(){
    return end_it->second;
}

Coord PeakReader::coordinates(){
    return Coord (end_it->second[0]->chr, end_it->second[0]->start, end_it->second[0]->end);
}

bool PeakReader::load(const std::string & full_filename){
    ifstream input_stream (full_filename);
    if (!input_stream) {
        cout << "Cannot open file " << full_filename << endl;
        return false;
    }
    string line;
    while (getline(input_stream, line)) {
        if (include_key(line, "#") || include_key(line, "start")) { // to filter commented lines and header
//            cout << "Filtered line: " << line << endl;
            continue;
        }
        PeakRecordPtr new_peak_record_ptr;
        try {
            PeakRecordPtr temp (new PeakRecord (line) ); // Could be redundant
            new_peak_record_ptr = temp;
        } catch (...){
//            cout << "Skipped line [" << line << "]" << endl;
            continue;
        }

//      [chromosome][peak start][peak end] = vector of PeakRecord
//      std::map<std::string, std::map<int, std::map<int, std::vector<PeakRecordPtr> > > > peak_data;
        std::string chr = new_peak_record_ptr->chr;
        int start = new_peak_record_ptr->start;
        int end = new_peak_record_ptr->end;
        if (! include_chr(chr)){
            std::vector <PeakRecordPtr> peak_ptr_vectror;
            std::map<int, std::vector<PeakRecordPtr> > end_layer_map;
            std::map<int, std::map<int, std::vector<PeakRecordPtr> > > start_layer_map;
            peak_ptr_vectror.push_back(new_peak_record_ptr);
            end_layer_map[end] = peak_ptr_vectror;
            start_layer_map[start] = end_layer_map;
            peak_data[chr] = start_layer_map;
        } else if (! include_start(chr, start)){
            std::vector <PeakRecordPtr> peak_ptr_vectror;
            std::map<int, std::vector<PeakRecordPtr> > end_layer_map;
            peak_ptr_vectror.push_back(new_peak_record_ptr);
            end_layer_map[end] = peak_ptr_vectror;
            peak_data[chr][start] = end_layer_map;
        } else if (! include_end(chr, start, end)){
            std::vector <PeakRecordPtr> peak_ptr_vectror;
            peak_ptr_vectror.push_back(new_peak_record_ptr);
            peak_data[chr][start][end] = peak_ptr_vectror;
        } else {
            peak_data[chr][start][end].push_back(new_peak_record_ptr);
        }
    }
    if (peak_data.size() == 0){
        // cout << "Couldn't load any data from file" < endl;
        return false;
    }
    reset();
    filename = full_filename;
    return true;
}

bool PeakReader::update (std::string chr, int start, int end, const GeneInfo & new_gene_info){
    if (include_end(chr, start, end)){
        for (int i = 0; i < peak_data[chr][start][end].size(); i++){
            peak_data[chr][start][end][i]->gene_info = new_gene_info;
        }
        return true;
    }
    return false;
}

void PeakReader::print(ostream& output_stream){
    reset();
    output_stream // header line
            << "refseq_id" << "\t"
            << "gene_id" << "\t"
            << "txStart" << "\t"
            << "txEnd" << "\t"
            << "strand" << "\t"
            << "chrom" << "\t"
            << "start" << "\t"
            << "end" << "\t"
            << "length" << "\t"
            << "abssummit" << "\t"
            << "pileup" << "\t"
            << "log10p" << "\t"
            << "foldenrich" << "\t"
            << "log10q" << "\t"
            << "region" << endl;
    do {
        for (int i = 0; i < value().size(); i++){
            output_stream << value()[i]->gene_info.refseq_id << "\t";
            output_stream << value()[i]->gene_info.gene_id << "\t";
            output_stream << value()[i]->gene_info.txStart << "\t";
            output_stream << value()[i]->gene_info.txEnd << "\t";
            output_stream << value()[i]->gene_info.strand << "\t";
            output_stream << value()[i]->chr << "\t";
            output_stream << value()[i]->start << "\t";
            output_stream << value()[i]->end << "\t";
            output_stream << value()[i]->length << "\t";
            output_stream << value()[i]->abs_summit << "\t";
            output_stream << value()[i]->pileup << "\t";
            output_stream << value()[i]->log10p << "\t";
            output_stream << value()[i]->fold_enrichment << "\t";
            output_stream << value()[i]->log10q << "\t";
            output_stream << value()[i]->gene_info.region << endl;
        };
    } while (next());
}

std::string PeakReader::get_filename(){
    return filename;
}

bool PeakReader::save (const std::string & output_filename){
    ofstream output_stream (output_filename);
    if (output_stream.is_open())
    {
        print(output_stream);
        output_stream.close();
        return true;
    }
    return false;
}

void PeakRecord::print(){
    cout << "PeakRecord" << endl;
    cout << "   start: " << start << endl;
    cout << "   end: " << end << endl;
    cout << "   length: " << length << endl;
    cout << "   abs_summit: " << abs_summit << endl;
    cout << "   pileup: " << pileup << endl;
    cout << "   log10p: " << log10p << endl;
    cout << "   fold_enrichment: " << fold_enrichment << endl;
    cout << "   log10q: " << log10q << endl;
    cout << "   name: " << name << endl;
    gene_info.print();
    cout << endl;
}

void GeneInfo::print(){
    cout << "GeneInfo" << endl;
    cout << "   refseq_id: " << refseq_id << endl;
    cout << "   gene_id: " << gene_id << endl;
    cout << "   txStart: " << txStart << endl;
    cout << "   txEnd: " << txEnd << endl;
    cout << "   strand: " << strand << endl;
    cout << "   region: " << region << endl;
    cout << endl;
}
