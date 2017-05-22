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
    if (not str_to_long(start, line_splitted[1])){
        throw ("PeakRecord constructor error");
    }
    // end
    if (not str_to_long(end, line_splitted[2])){
        throw ("PeakRecord constructor error");
    }
    // length
    if (not str_to_int(length, line_splitted[3])){
        throw ("PeakRecord constructor error");
    }
    // abs_summit
    if (not str_to_long(abs_summit, line_splitted[4])){
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
    load (full_filename);
}


bool PeakReader::include_chr(const std::string chr){
    std::map<std::string, std::map<long, std::map<long, std::vector<PeakRecordPtr> > > >::iterator it = peak_data.find(chr);
    if(it != peak_data.end()) {
        return true;
    }
    return false;
}

bool PeakReader::include_start(const std::string chr, long start){
    if (include_chr(chr)){
        std::map<long, std::map<long, std::vector<PeakRecordPtr> > >::iterator it = peak_data[chr].find(start);
        if(it != peak_data[chr].end()) {
            return true;
        }
    }
    return false;
}

bool PeakReader::include_end(const std::string chr, long start, long end){
    if ( include_start(chr,start) ){
        std::map<long, std::vector<PeakRecordPtr> >::iterator it = peak_data[chr][start].find(end);
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
    std::map<std::string, std::map<long, std::map<long, std::vector<PeakRecordPtr> > > >::iterator chr_it_copy = chr_it;
    std::map<long, std::map<long, std::vector<PeakRecordPtr> > > start_it_copy = start_it;
    std::map<long, std::vector<PeakRecordPtr> > end_it_copy = end_it;
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

PeakRecordPtr PeakReader::value(){
    return end_it->second;
}

Coord PeakReader::coordinates(){
    return Coord (end_it->second[0].chr, end_it->second[0].start, end_it->second[0].end);
}

PeakReader::load(const std::string & full_filename){
    ifstream input_stream (full_filename);
    if (!input_stream) {
        cout << "Cannot open file " << full_filename << endl;
        return false;
    }
    string line;
    while (getline(input_stream, line)) {
        if (include_key(line, "#") || include_key(line, "start")) { // to filter commented lines and header
            cout << "Filtered line: " << line << endl;
            continue;
        }
        PeakRecordPtr new_peak_record_ptr;
        try {
            PeakRecordPtr temp (new PeakRecord (line) ); // Could be redundant
            new_peak_record_ptr = temp;
        } catch (...){
            cout << "Skipped line [" << line << "]" << endl;
            continue;
        }

//      [chromosome][peak start][peak end] = vector of PeakRecord
//      std::map<std::string, std::map<long, std::map<long, std::vector<PeakRecordPtr> > > > peak_data;
        std::string chr = new_peak_record.chr;
        std::string start = new_peak_record.start;
        std::string end = new_peak_record.end;
        if (! include_chr(chr)){
            std::vector <PeakRecordPtr> peak_ptr_vectror;
            std::map<long, std::vector<PeakRecordPtr> > end_layer_map;
            std::map<long, std::map<long, std::vector<PeakRecordPtr> > > start_layer_map;
            peak_ptr_vectror.push_back(new_peak_record_ptr);
            end_layer_map[end] = peak_ptr_vectror;
            start_layer_map[start] = end_layer_map;
            peak_data[chr] = start_layer_map;
        } else if (! include_start(chr, start)){
            std::vector <PeakRecordPtr> peak_ptr_vectror;
            std::map<long, std::vector<PeakRecordPtr> > end_layer_map;
            eak_ptr_vectror.push_back(new_peak_record_ptr);
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
    reset();
}

bool PeakReader::update (std::string chr, long start, long end, GeneInfo new_gene_info){
    if (include_end(chr, start, end)){
        for (int i = 0; peak_data[chr][start][end].length(); i++){
            peak_data[chr][start][end][i].gene_info = new_gene_info;
        }
        return true;
    }
    return false;
}
