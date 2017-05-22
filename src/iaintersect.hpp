#ifndef _iaintersect_
#define _iaintersect_


#include <config.hpp>

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
        qint64 txStart;
        qint64 txEnd;
        Annotation(){};
        Annotation(QString c,
                   QChar s,
                   int ec,
                   QStringList eS,
                   QStringList eE,
                   QString n,
                   QString n2,
                   qint64 txS,
                   qint64 txE):
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

typedef qint64 t_genome_coordinates;
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
