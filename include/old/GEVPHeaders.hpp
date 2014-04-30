/*
 * GEVPHeaders.hpp
 *
 *  Created on: Jun 14, 2013
 *      Author: Thibaut Metivet
 */

#ifndef GEVP_HEADERS_HPP_
#define GEVP_HEADERS_HPP_

#include "LQCDA.hpp"

    class GEVPCorrelatorHeader
    {
    public:
	struct Smearing {
	    enum class Smear { S, L, U };
	    Smear Sm;
	    Smearing() : Sm(Smear::U) {}
	    Smearing(Smear sm) : Sm(sm) {}
	    Smearing(const std::string& s) {
		if(s == "S") Sm = Smear::S;
		if(s == "L") Sm = Smear::L;
	    }
	};
	struct Interpolator {
	    enum class Interpol { Vi, P, U };
	    Interpol I;
	    Interpolator() : I(Interpol::U) {}
	    Interpolator(Interpol i) : I(i) {}
	    Interpolator(const std::string& s) {
		if(s == "Vi") I = Interpol::Vi;
		if(s == "P") I = Interpol::P;
	    }
	};  

	Smearing SourceSmearing, SinkSmearing;
	Interpolator Source, Sink;
	unsigned int Peq;
	int Pola;

    public:
	GEVPCorrelatorHeader() {}
	GEVPCorrelatorHeader(const char* srcsm, const char* snksm,
			     const char* src, const char* snk,
			     unsigned int peq, int pola) :
	    SourceSmearing(std::string(srcsm)), SinkSmearing(std::string(snksm)),
	    Source(std::string(src)), Sink(std::string(snk)),
	    Peq(peq), Pola(pola) {}

	GEVPCorrelatorHeader(const std::string& s);
    
    };

    bool operator== (const GEVPCorrelatorHeader& h1, const GEVPCorrelatorHeader& h2);

#endif	// GEVP_HEADERS_HPP_
