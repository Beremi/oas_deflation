#ifndef _VTK_EXPORTER_H
#define _VTK_EXPORTER_H

#include "data_exporter.h"


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// VTK EXPORTERS
class VTKExporter : public DataExporter
{
private:

public:
    VTKExporter(unsigned dimension) : DataExporter(dimension) {};
    virtual ~VTKExporter() {};
    virtual void readFromLine(istringstream &iss);
    virtual void giveFileName(unsigned step, char *buffer) const;
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions) const = 0;
protected:
    unsigned cell_data_size;
};


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT ELEMENTS TO VTK (VTU)
class VTKElementExporter : public VTKExporter
{
private:
    ElementContainer *elems;
    NodeContainer *nodes;
public:
    VTKElementExporter(ElementContainer *e, NodeContainer *n, unsigned dimension) : VTKExporter(dimension) { elems = e; nodes = n; };
    ~VTKElementExporter() {};
    // virtual void readFromLine(istringstream &iss);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT RIGID BODIES TO VTK (VTU)
// WORKS ONLY FOR RIGID BODIES IN 2D
class VTKRB2DExporter : public VTKExporter
{
private:
    ElementContainer *elems;
    NodeContainer *nodes;
public:
    VTKRB2DExporter(ElementContainer *e, NodeContainer *n, unsigned dimension) : VTKExporter(dimension) { elems = e; nodes = n; };
    ~VTKRB2DExporter() {};
    // virtual void readFromLine(istringstream &iss, unsigned dimension);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions) const;
protected:
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
// EXPORT RIGID contacts TO VTK (VTU)
// TODO add posibility to export only certain nodes
class VTKRCExporter : public VTKExporter
{
private:
    ElementContainer *elems;
    NodeContainer *nodes;
public:
    VTKRCExporter(ElementContainer *e, NodeContainer *n, unsigned dimension) : VTKExporter(dimension) { elems = e; nodes = n; };
    ~VTKRCExporter() {};
    // virtual void readFromLine(istringstream &iss, unsigned dimension);
    virtual void exportData(unsigned step, const Vector &DoFs, const Vector &reactions) const;
protected:
};

#endif /* _EXPORTER_C_H */
