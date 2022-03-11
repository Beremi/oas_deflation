#ifndef _EXPORTER_M_H
#define _EXPORTER_M_H

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
    virtual void giveFileName(unsigned step, char *buffer) const;
    virtual void exportData(unsigned step, const MyVector &DoFs, const MyVector &reactions, fs :: path resultDir) const;
protected:
};



#endif /* _EXPORTER_M_H */
