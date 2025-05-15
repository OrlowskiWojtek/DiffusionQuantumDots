#ifndef STRUCT_WALKER_HPP
#define STRUCT_WALKER_HPP

#include <array>
#include <vector>

struct walker {
    std::array<double, 3> cords;

    const double &x() const { return cords[0]; };
    const double &y() const { return cords[1]; };
    const double &z() const { return cords[2]; };

    walker() = default;
    walker(double xval, double yval, double zval){
        cords[0] = xval;
        cords[1] = yval;
        cords[2] = zval;
    }
};

// electron_walker representing simple walker in many dimensions systems
// it is made in order to increase number of dimensions without 
// increasing number of dimensions in walker_struct
// for 1 electron_walker vector size is 1
using electron_walker = std::vector<walker>;

#endif // ndef STRUCT_WALKER_HPP
