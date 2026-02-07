#ifndef _EXPORTER_MODEL_H
#define _EXPORTER_MODEL_H

#include "data_exporter.h"

class ElementStatsExporter : public DataExporter
{
private:
    bool remove_previous;
    std :: string last_saved_file;
    std :: vector< unsigned >elems_to_save;
    ElementContainer *elems;
public:
    ElementStatsExporter(ElementContainer *e, unsigned dimension) : DataExporter(dimension) { elems = e; };
    virtual ~ElementStatsExporter() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual std :: string giveFileName(unsigned step, int iteration) const override;
    virtual void exportData(unsigned step, int iteration, fs :: path resultDir) const;
protected:
};



#endif /* _EXPORTER_MODEL_H */
