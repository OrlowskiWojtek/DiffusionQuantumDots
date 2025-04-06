#ifndef STRUCT_WALKER_HPP
#define STRUCT_WALKER_HPP

#include <array>

struct walker {
    std::array<double, 3> cords;

    double &x() { return cords[0]; };
    double &y() { return cords[1]; };
    double &z() { return cords[2]; };

    walker() = default;
    walker(double xval, double yval, double zval){
        cords[0] = xval;
        cords[1] = yval;
        cords[2] = zval;
    }
};

#endif // ndef STRUCT_WALKER_HPP
