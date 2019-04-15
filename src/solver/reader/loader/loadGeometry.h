/**
 * @Author: jose
 * @Date:   2019-04-05T19:01:08+02:00
 * @Last modified by:   jose
 * @Last modified time: 2019-04-10T19:18:25+02:00
 */

// #include <iostream>
#include <stdio.h>
#include <vector>
#include "../../common/point.h"
#include "../../prepro/nodes/node.h"
#include "../../prepro/nodes/vertex.h"
#include "../../prepro/elements/element_master.h"
#include "../../common/superintendent.h"

#ifndef _LOADGEOMETRY_H
#define _LOADGEOMETRY_H

void loadInput(Superintendent &SUP, std::vector<std::vector<Element *> > &ELEMENTS, std::vector<Node> &Nodes, std::vector<Vertex> &Vertices);


#endif	/* _LOADGEOMETRY_H */
