#ifndef UNIT_HANDLER_HPP
#define UNIT_HANDLER_HPP

// Class for easy conversion between SI and atomic units
class UnitHandler{
    // DEFAULT is defined as SI, but nm instead of m, and meV instead of eV
    static constexpr double bohr_radius_nm  = 0.0529;
    static constexpr double hartree_meV     = 27211.6;

public:
    enum ConvMode{ TO_AU, TO_DEFAULT};

    static constexpr double length(ConvMode conv, double value){
        switch(conv){
            case TO_AU: return value / bohr_radius_nm; break;
            case TO_DEFAULT: return value * bohr_radius_nm; break;
            default: return 0; break;
        }
    }
    
    static constexpr double energy(ConvMode conv, double value){
        switch(conv){
            case TO_AU: return value / hartree_meV; break;
            case TO_DEFAULT: return value * hartree_meV; break;
            default: return 0; break;
        }
    }
};

#endif
