#include "iaintersect.hpp"

IAIntersect::IAIntersect(QObject *parent):
    QObject(parent)
{
    promoter=gArgs().getArgs("promoter").toInt();
    upstream=gArgs().getArgs("upstream").toInt();
}

IAIntersect::~IAIntersect()
{
}

AnnotationMap& AnnotationMap::operator+=(const AnnotationMap& other)
{
    for(auto key : other.data.keys()) {
        this->data.insertMulti(key, other.data.value(key));
    }
    return *this;
}

bool AnnotationMap::operator==(const AnnotationMap& other) const
{
    return this->data == other.data;
}


void IAIntersect::fillUpAnnotation( ) {

    if(gArgs().getArgs("annotation").toString().isEmpty()) {
        throw "Set gene annotaion file";
    }

    string annotation_filename = gArgs().getArgs("annotation").toString().toStdString();
    ifstream annotation_stream (annotation_filename.c_str());
    if (!annotation_stream) {
        throw "Failed to open annotaion file";
    }

    std::string line;
    getline(annotation_stream, line);
    std::map<std::string,short> column_map;
    if (line.length() > 0 && line.at(0) == '#'){
        QStringList header_splitted = QString::fromStdString(line).split(QChar('\t'));
        for (int i = 0; i < header_splitted.length(); i++){
            column_map[header_splitted.at(i).toStdString()] = i;
        }
    } else {
        column_map["exonCount"] = 8;
        column_map["exonStarts"] = 9;
        column_map["exonEnds"] = 10;
        column_map["name"] = 1;
        column_map["name2"] = 12;
        column_map["chrom"] = 2;
        column_map["strand"] = 3;
        column_map["txStart"] = 4;
        column_map["txEnd"] = 5;
    }

    while(getline(annotation_stream, line)) {
        QStringList line_splitted = QString::fromStdString(line).split(QChar('\t'));

        QStringList q_starts = line_splitted.at(column_map["exonStarts"]).split(QChar(','));
        QStringList q_ends = line_splitted.at(column_map["exonEnds"]).split(QChar(','));
        int exCount = line_splitted.at(column_map["exonCount"]).toInt();
        QString chr = line_splitted.at(column_map["chrom"]);
        QChar strand = line_splitted.at(column_map["strand"]).at(0);
        QString iso_name = line_splitted.at(column_map["name"]);
        QString g_name = line_splitted.at(column_map["name2"]);
        quint64 txStart = line_splitted.at(column_map["txStart"]).toInt()+1;
        quint64 txEnd = line_splitted.at(column_map["txEnd"]).toInt();



        if(gArgs().getArgs("sam_ignorechr").toString().contains(chr)) {
            continue;
        }

        annotationPtr p( new Annotation(chr,strand,exCount,q_starts,q_ends,iso_name,g_name,txStart,txEnd) );

        AnnotationMap aptr;
        aptr.data.insert(coord_key(txStart,txEnd), p);

        if(strand=='+') {
            if(promoter+upstream >txStart)
                txStart=1;
            else
                txStart-=(promoter+upstream);
        } else {
            txEnd+=(promoter+upstream);
        }
        annotation[chr].add(make_pair(bicl::discrete_interval<t_genome_coordinates>::closed(txStart,txEnd),aptr));
        for(int j=0;j<exCount;j++) {
            int s=q_starts.at(j).toInt(),e=q_ends.at(j).toInt();
            annotation[chr].add(make_pair(bicl::discrete_interval<t_genome_coordinates>::closed(s,e),aptr));
        }

    }

}


void IAIntersect::start() {
    this->fillUpAnnotation();

    if(gArgs().getArgs("in").toString().isEmpty()){
        throw "Set input file";
    }
    if(gArgs().getArgs("out").toString().isEmpty()){
        throw "Set output file";
    }

    PeakReader peak_reader (gArgs().getArgs("in").toString().toStdString());

    peak_reader.reset();
    do {
//        cout << peak_reader.coordinates().chr << '\t' << peak_reader.coordinates().start << '\t' << peak_reader.coordinates().end << endl;
        QStringList GENE_NAME;
        QStringList REFSEQ_NAME;
        int txStart;
        int txEnd;
        QChar strand;

        QString chr=QString::fromStdString(peak_reader.coordinates().chr);
        int start=peak_reader.coordinates().start;
        int end=peak_reader.coordinates().end;

        if(!annotation.keys().contains(chr)) continue;
        bicl::discrete_interval<t_genome_coordinates> current_region=bicl::discrete_interval<t_genome_coordinates>::closed(start,end);
        pair<chrom_coverage::iterator,chrom_coverage::iterator> pi=annotation[chr].equal_range(current_region);

        if(pi.first==pi.second) {
            chrom_coverage::iterator itl=annotation[chr].lower_bound(current_region);
            chrom_coverage::iterator itu=annotation[chr].upper_bound(current_region);
            chrom_coverage::iterator it;
            if(itl==annotation[chr].end())
                itl--;
            if(itu==annotation[chr].end())
                itu--;
            if(itl==annotation[chr].end() && itu==annotation[chr].end())
                continue;
            if(itl!=annotation[chr].end() && itu!=annotation[chr].end()) {
                if(itl->first.lower()-end < start-itu->first.upper())
                    it=itl;
                else
                    it=itu;
            }
            else if(itu==annotation[chr].end())
                it=itl;
            else
                it=itu;

            QMapIterator<coord_key, annotationPtr> i((*it).second.data);

            while(i.hasNext()){
                annotationPtr p(i.next().value());
                if(GENE_NAME.isEmpty()) {
                    txStart=p->txStart;
                    txEnd=p->txEnd;
                    strand=p->strand;
                }
                if(!GENE_NAME.contains(p->name2))
                    GENE_NAME<<p->name2;
                if(!REFSEQ_NAME.contains(p->name))
                    REFSEQ_NAME<<p->name;
            }
            GENE_NAME.sort();
            REFSEQ_NAME.sort();
            QString gene_name,refseq_name;
            gene_name=GENE_NAME.value(0);
            refseq_name=REFSEQ_NAME.value(0);
            for(int i=1;i<GENE_NAME.count();i++){
                gene_name+=","+GENE_NAME.value(i);
                refseq_name+=","+REFSEQ_NAME.value(i);
            }
            if(!peak_reader.update(chr.toStdString(), start, end, GeneInfo(refseq_name.toStdString(), gene_name.toStdString(), txStart, txEnd, strand.toLatin1(), "intergenic")) ) {
                throw "Error updating islands";
            }
            continue;
        }
        QMap<int,AnnotationMap> result;

        while(pi.first!=pi.second){
            //chrom_coverage::interval_type itv = bicl::key_value<chrom_coverage >(pi.first);
//            qDebug()<<chr<<"["<<itv.lower()<<":"<<itv.upper()<<"]=["<<start<<":"<<end<<"]";

            QMapIterator<coord_key, annotationPtr> i((*(pi.first)).second.data);
            while(i.hasNext()){
                annotationPtr p(i.next().value());
                bicl::discrete_interval<t_genome_coordinates> _promoter;
                bicl::discrete_interval<t_genome_coordinates> _upstream;
                bicl::discrete_interval<t_genome_coordinates> _gene;

                if(p->strand=='+') {
                    int beg=p->txStart-promoter;
                    _promoter =bicl::discrete_interval<t_genome_coordinates>::closed( ( beg< 0 ) ? 0: (beg),p->txStart+promoter);
                    beg-=upstream;
                    _upstream=bicl::discrete_interval<t_genome_coordinates>::closed( (beg)? 0: (beg),p->txStart-promoter);
                } else {
                    int beg=p->txEnd-promoter;
                    _promoter =bicl::discrete_interval<t_genome_coordinates>::closed( (beg <0 )?0:(beg),p->txEnd+promoter);
                    _upstream =bicl::discrete_interval<t_genome_coordinates>::closed( p->txEnd+promoter,p->txEnd+promoter+upstream);
                }
                _gene =bicl::discrete_interval<t_genome_coordinates>::closed(p->txStart,p->txEnd);


                //is promoter
                if(bicl::intersects(current_region,_promoter)) {
                    result[1].data.insert(coord_key(p->txStart,p->txEnd), p);
                    continue;
                }

                if(result.contains(1))//skeep all others if promoter exists
                    continue;
                //is intron
                if(bicl::intersects(current_region,_gene)) {
                    result[3].data.insert(coord_key(p->txStart,p->txEnd), p);
                }

                //is Exon
                bicl::interval_set<t_genome_coordinates> _exons;
                for(int j=0;j<p->exonCount;j++) {
                    quint64 s=p->exonStarts.at(j).toInt(),e=p->exonEnds.at(j).toInt();
                    _exons.add(bicl::discrete_interval<t_genome_coordinates>::closed(s,e));
                }

                if(bicl::intersects(_exons,current_region)) {
                    result[2].data.insert(coord_key(p->txStart,p->txEnd), p);
                    continue;
                }
                if(result.contains(2) || result.contains(3))//skeep all others
                    continue;

                //is upstream
                if(bicl::intersects(current_region,_upstream)) {
                    result[4].data.insert(coord_key(p->txStart,p->txEnd), p);
                    //                    qDebug()<<chr<<"["<<_promoter.lower()<<":"<<_promoter.upper()<<"]=["<<current_region.lower()<<":"<<current_region.upper()<<"]";
                    continue;
                }
            }//while i.hasnext
            pi.first++;
        }//while end!=begin

        //update db
        QMap<int,AnnotationMap >::Iterator mit=result.upperBound(0);
        QString region;
        switch(mit.key()) {
            case 1:
                region="promoter";
                break;
            case 2:
                region="exon";
                break;
            case 3:
                region="intron";
                break;
            case 4:
                region="upstream";
                break;
        }

        QMapIterator<coord_key, annotationPtr> i(mit.value().data);
        while(i.hasNext()){
            annotationPtr p(i.next().value());
            if(GENE_NAME.isEmpty()) {
                txStart=p->txStart;
                txEnd=p->txEnd;
                strand=p->strand;
            }
            if(!GENE_NAME.contains(p->name2))
                GENE_NAME<<p->name2;
            if(!REFSEQ_NAME.contains(p->name))
                REFSEQ_NAME<<p->name;
        }
        GENE_NAME.sort();
        REFSEQ_NAME.sort();
        QString gene_name,refseq_name;
        gene_name=GENE_NAME.value(0);
        refseq_name=REFSEQ_NAME.value(0);
        for(int i=1;i<GENE_NAME.count();i++){
            gene_name+=","+GENE_NAME.value(i);
            refseq_name+=","+REFSEQ_NAME.value(i);
        }
        if(!peak_reader.update(chr.toStdString(), start, end, GeneInfo(refseq_name.toStdString(), gene_name.toStdString(), txStart, txEnd, strand.toLatin1(), region.toStdString())) ){
            throw "Error updating islands";
        }
    } while (peak_reader.next());

    if ( !peak_reader.save( gArgs().getArgs("out").toString().toStdString() ) ){
        throw "Error export results";
    }

    emit finished();
}
