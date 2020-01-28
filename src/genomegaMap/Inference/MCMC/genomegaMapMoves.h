/*  Copyright 2018 Daniel Wilson.
 *
 *  Part of the genomegaMap library.
 *
 *  The genomegaMap library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The genomegaMap library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the genomegaMap library. If not, see <http://www.gnu.org/licenses/>.
 */
/*
 *  genomegaMapMoves.h
 *  gcat
 *
 *  Created by Daniel Wilson on 29/04/2010.
 *
 */
#ifndef _GENOMEGAMAP_MCMC_MOVES_H_
#define _GENOMEGAMAP_MCMC_MOVES_H_
#include <Inference/MCMC/MCMC.h>

using namespace gcat;

namespace genomegaMap {
	
class Codon61SequenceGibbsSampler : public MCMC_move {
public:
	// Constructor
	Codon61SequenceGibbsSampler(MCMC* mcmc, vector< std::string > &target, const double weight, std::string type="Codon61SequenceGibbsSampler_move");
	// Go!!!
	void go();
};
	
} // namespace genomegaMap


#endif // _GENOMEGAMAP_MCMC_MOVES_H_
