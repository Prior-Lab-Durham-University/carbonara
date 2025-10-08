#include "point.h"

point::point(double x, double y, double z) : X(x), Y(y), Z(z), norm(0.0) {}
point::point() : X(0.0), Y(0.0), Z(0.0), norm(0.0) {}

point::point(const std::string& triplet) : X(0.0), Y(0.0), Z(0.0), norm(0.0) {
    std::stringstream ss(triplet);
    double d = 0.0;
    std::vector<double> dv;
    while (ss >> d) dv.push_back(d);
    if (dv.size() >= 3) { X = dv[0]; Y = dv[1]; Z = dv[2]; }
}

double point::getX() const { return X; }
double point::getY() const { return Y; }
double point::getZ() const { return Z; }

void point::setX(double xv) { X = xv; }
void point::setY(double yv) { Y = yv; }
void point::setZ(double zv) { Z = zv; }

double point::length() const { return std::sqrt(X*X + Y*Y + Z*Z); }

point point::sum(const point& b) const { return point(b.getX() + X, b.getY() + Y, b.getZ() + Z); }
point point::dif(const point& b) const { return point(X - b.getX(), Y - b.getY(), Z - b.getZ()); }

point point::cross(const point& p2) const {
    return point(Y * p2.getZ() - Z * p2.getY(),
                 Z * p2.getX() - X * p2.getZ(),
                 X * p2.getY() - Y * p2.getX());
}

double point::dotprod(const point& p2) const { return X*p2.getX() + Y*p2.getY() + Z*p2.getZ(); }

double point::scalarTriple(const point& p1, const point& p2, const point& p3) {
    point cp = p2.cross(p3);
    return p1.dotprod(cp);
}

void point::scalarMult(double a) {
    X *= a; Y *= a; Z *= a;
}

bool point::checkEqual(const point& p2) const {
    return (X == p2.getX() && Y == p2.getY() && Z == p2.getZ());
}

void point::normalise() {
    norm = std::sqrt(X*X + Y*Y + Z*Z);
    if (norm != 0.0) { X /= norm; Y /= norm; Z /= norm; }
}

void point::znormalise() {
    if (std::abs(Z) > 1e-8) {
        X /= Z; Y /= Z; Z = 1.0;
    }
}

double point::pairDist(const point& p) const {
    // legacy behavior: returns 0.0 if almost equal, 1.0 otherwise
    double xd = X - p.getX();
    double yd = Y - p.getY();
    double zd = Z - p.getZ();
    if (std::abs(xd) < 1e-8 && std::abs(yd) < 1e-8 && std::abs(zd) < 1e-8) {
        return 0.0;
    } else {
        return 1.0;
    }
}

bool point::isNonzero() const {
    return !(std::abs(X) < 1e-8 && std::abs(Y) < 1e-8 && std::abs(Z) < 1e-8);
}

void point::printPoint() const { std::cout << X << " " << Y << " " << Z << "\n"; }

point point::operator+(const point& p) const { return point(X + p.getX(), Y + p.getY(), Z + p.getZ()); }
point point::operator-(const point& p) const { return point(X - p.getX(), Y - p.getY(), Z - p.getZ()); }

point point::operator*(double d) const { return point(d*X, d*Y, d*Z); }
point point::operator/(double d) const { return point(X/d, Y/d, Z/d); }

double point::eDist(const point& p2) const {
    double xd = X - p2.getX();
    double yd = Y - p2.getY();
    double zd = Z - p2.getZ();
    return std::sqrt(xd*xd + yd*yd + zd*zd);
}
