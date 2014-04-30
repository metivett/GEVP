/*
 * GEVP.cpp
 *
 *  Created on: Jun 14, 2013
 *      Author: Thibaut Metivet
 */

#include "GEVPHeaders.hpp"
#include "boost/regex.hpp"

using namespace boost;

GEVPCorrelatorHeader::GEVPCorrelatorHeader(const std::string& s)
{
    smatch sm;
    regex rx("#\\s+(S|L)\\s+(S|L)\\s+(Vi|P)(Vi|P)_Peq(\\d)(pola[0-2])?\\s+0\\s+0");
    regex_match(s, sm, rx);
    SourceSmearing = Smearing(sm[1]);
    SinkSmearing = Smearing(sm[2]);
    Source = Interpolator(sm[3]);
    Sink = Interpolator(sm[4]);
    Peq = atoi(sm[5].str().c_str());
    Pola = (!sm[6].str().empty()) ? (int)(sm[6].str()[4] - 48) : (-1);
}

bool operator== (const GEVPCorrelatorHeader::Smearing& s1, const GEVPCorrelatorHeader::Smearing& s2)
{
    return s1.Sm == s2.Sm;
}

bool operator== (const GEVPCorrelatorHeader::Interpolator& i1, const GEVPCorrelatorHeader::Interpolator& i2)
{
    return i1.I == i2.I;
}

bool operator== (const GEVPCorrelatorHeader& h1, const GEVPCorrelatorHeader& h2)
{
    return h1.SourceSmearing == h2.SourceSmearing
	&& h1.SinkSmearing == h2.SinkSmearing
	&& h1.Source == h2.Source
	&& h1.Sink == h2.Sink
	&& h1.Peq == h2.Peq;
}
