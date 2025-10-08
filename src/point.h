#ifndef PNT_H
#define PNT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>

class point {
public:
    point(double x, double y, double z);
    point();
    explicit point(const std::string& triplet);

    // getters
    double getX() const;
    double getY() const;
    double getZ() const;

    // setters (mutating)
    void setX(double xv);
    void setY(double yv);
    void setZ(double zv);

    // geometry
    double length() const;
    point  sum(const point& b) const;
    point  dif(const point& b) const;
    point  cross(const point& p2) const;
    double dotprod(const point& p2) const;
    static double scalarTriple(const point& p1, const point& p2, const point& p3);

    // comparisons / checks
    bool   checkEqual(const point& p2) const;
    bool   isNonzero() const;

    // mutating ops
    void   scalarMult(double a);
    void   normalise();
    void   znormalise();

    // misc
    void   printPoint() const;

    // operators (non-mutating)
    point  operator+(const point& p) const;
    point  operator-(const point& p) const;
    point  operator*(double d) const;
    point  operator/(double d) const;

    // distances
    double eDist(const point& p2) const;

    // NOTE: legacy: returns 0.0/1.0 like before (true => 1.0)
    double pairDist(const point& p) const;

private:
    double X, Y, Z, norm;
};

#endif
