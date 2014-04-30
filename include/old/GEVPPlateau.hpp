/*
 * GEVPPlateau.hpp
 *
 *  Created on: Jun 18, 2013
 *      Author: Thibaut Metivet
 */

#ifndef GEVP_PLATEAU_HPP_
#define GEVP_PLATEAU_HPP_

#include <iostream>

    struct GEVPPlateau {
	double E1, E2;
	GEVPPlateau() {}
	GEVPPlateau(double e1, double e2) :
	    E1(e1), E2(e2) {}
    };


    inline std::ostream& operator<< (std::ostream& out, const GEVPPlateau& plat) {
	out << "(E1 = " << plat.E1 << ", E2 = " << plat.E2 << ")" << std::endl;
	return out;
    }

    inline GEVPPlateau operator+ (const GEVPPlateau& p1, const GEVPPlateau& p2) {
	return GEVPPlateau(p1.E1 + p2.E1, p1.E2 + p2.E2);
    }

    inline GEVPPlateau operator- (const GEVPPlateau& p1, const GEVPPlateau& p2) {
	return GEVPPlateau(p1.E1 - p2.E1, p1.E2 - p2.E2);
    }

    inline GEVPPlateau operator* (const GEVPPlateau& p1, const GEVPPlateau& p2) {
	return GEVPPlateau(p1.E1 * p2.E1, p1.E2 * p2.E2);
    }

    inline GEVPPlateau operator+ (const GEVPPlateau& p1, double d) {
	return GEVPPlateau(p1.E1 + d, p1.E2 + d);
    }

    inline GEVPPlateau operator- (const GEVPPlateau& p1, double d) {
	return GEVPPlateau(p1.E1 - d, p1.E2 - d);
    }

    inline GEVPPlateau operator* (const GEVPPlateau& p, double d) {
	return GEVPPlateau(p.E1 * d, p.E2 * d);
    }

    inline GEVPPlateau operator/ (const GEVPPlateau& p, double d) {
	return GEVPPlateau(p.E1 / d, p.E2 / d);
    }

inline GEVPPlateau sqrt(const GEVPPlateau& p) {
    return GEVPPlateau(std::sqrt(p.E1), std::sqrt(p.E2));
}


#endif	// GEVP_PLATEAU_HPP_
