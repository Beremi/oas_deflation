/**
 * @Author: jose
 * @Date:   2019-04-10T15:48:44+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T17:36:28+02:00
 */

#include "../../common/point.h"
#include <vector>

#ifndef _NODE_H
#define _NODE_H


class Node {
private:
    int id;  // NOTE is it necessary? or only order in vector is enough?
    Point N;  // center of voronoi cell
    double r;  // radius in case of power tessellation
public:
    Node();
    Node(const std::string &str);  // construct node form line in input file nodes txt
    Node(std::istringstream &iss);
    virtual ~Node();

    int get_id() const;
    void set_id(const int &i);

    Point get_node() const;
    void set_node(const Point &P);

    double get_radius() const;
    void set_radius(const double &a);


protected:

};

void readNodes(const std::string &inputfile, std::vector<Node> &Nodes);

#endif /* _NODE_H */
