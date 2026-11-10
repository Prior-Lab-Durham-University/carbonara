#ifndef SKMT_H
#define SKMT_H

#include "point.h"
#include <iostream>
#include <vector>

class skmt{
    public:
        skmt();
        template <typename T> int sgn(T val);
        int signedTetraVolume(point&a,point&b,point&c,point&d);
        bool intersectLineTriangle(point&q1, point&q2, point&p1, point&p2, point&p3);
        std::vector<point> getSKMTCurve(ktlMolecule&molin);
};
#endif
