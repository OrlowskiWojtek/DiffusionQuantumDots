#ifndef ABSTRACT_ORBITAL_HPP
#define ABSTRACT_ORBITAL_HPP

class AbstractOrbital{
public:
    virtual void print() = 0;
    virtual ~AbstractOrbital() = default;

protected:
    virtual void print_test_to_file() = 0;
};

#endif
