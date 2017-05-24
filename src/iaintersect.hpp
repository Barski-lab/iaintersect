#ifndef _iaintersect_
#define _iaintersect_


#include <config.hpp>
#include "peak_reader.hpp"

#ifndef FSTM
#define FSTM IAIntersect
#endif

struct Annotation {
        QString chrom;
        QChar strand;
        int exonCount;
        QStringList exonStarts;
        QStringList exonEnds;
        QString name;
        QString name2;
        int txStart;
        int txEnd;
        Annotation(){};
        Annotation(QString c,
                   QChar s,
                   int ec,
                   QStringList eS,
                   QStringList eE,
                   QString n,
                   QString n2,
                   int txS,
                   int txE):
            chrom(c),
            strand(s),
            exonCount(ec),
            exonStarts(eS),
            exonEnds(eE),
            name(n),
            name2(n2),
            txStart(txS),
            txEnd(txE){};

};

typedef int t_genome_coordinates;
typedef QSharedPointer<Annotation> annotationPtr;
typedef bicl::split_interval_map<t_genome_coordinates,QSet< annotationPtr > > chrom_coverage;

class IAIntersect: public QObject
{
        Q_OBJECT
    private:
        quint64 promoter,upstream;
        QMap<QString, chrom_coverage > annotation;
        void getRecordInfo();
        void fillUpAnnotation();
    public slots:
        void start(void);

    signals:
        void finished(void);

    public:
        IAIntersect(QObject* parent=0);
        ~IAIntersect();
};

#endif
