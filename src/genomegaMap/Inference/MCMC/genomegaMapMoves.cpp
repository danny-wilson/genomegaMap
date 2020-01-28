/*  Copyright 2018 Daniel Wilson.
 *
 *  Part of the omegaMap library.
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
 *  genomegaMapMoves.cpp
 *  gcat
 *
 *  Created by Daniel Wilson on 29/04/2010.
 *
 */
#include <genomegaMap/Inference/MCMC/genomegaMapMoves.h>
#include <genomegaMap/RandomVariables/Codon61Sequence.h>

using namespace gcat;

namespace genomegaMap {
	
Codon61SequenceGibbsSampler::Codon61SequenceGibbsSampler(MCMC* mcmc, vector< string > &target, const double weight, string type) : MCMC_move(mcmc,target,weight,type) {
}

void Codon61SequenceGibbsSampler::go() {
	Codon61SequenceRV& rv = *((Codon61SequenceRV*)_target[0]);
	// Draw the site uniformly at random
	const int pos = _ran->discrete(0,rv.length()-1);
	mydouble old_likelihood = _mcmc->likelihood();
	// Cycle through all possible proposals, and accept them (except the last one).
	// Then draw the actual move with probability proportional to its joint likelihood. If not the last one, accept the last one and propose the one drawn.
	// (This has a slightly efficiency gain if it is the last one as it does not need to be re-proposed.) Finally, go through the usual motions.
	vector< mydouble > jlik(61);
	mydouble sumjlik = 0.0;
	int i;
	for(i=0;i<61;i++) {
		rv.change_value(pos,i,Variable::_PROPOSE);
		jlik[i] = _mcmc->update_likelihood();
		sumjlik += jlik[i];
		// Acceptance (rather than rejection) allows more efficient update algorithms to be implemented
		rv.change_value(pos,i,Variable::_ACCEPT);
	}
	// Draw the proposal
	double U = _ran->U();
	for(i=0;i<61;i++) {
		if((U -= (jlik[i]/sumjlik).todouble())<=0.0) break;
	}
	if(i==61) error("Codon61SequenceGibbsSampler::go(): did not draw move correctly");
	// For efficiency...
	if(i<60) {
		rv.change_value(pos,i,Variable::_PROPOSE);
		mydouble checklik = _mcmc->update_likelihood();
		if(abs(checklik.LOG()-jlik[i].LOG())>1e-6) {
			stringstream errMsg;
//			errMsg << "Codon61SequenceGibbsSampler::go(): inconsistent likelihood calculation. Was ";
//			errMsg << jlik[i].LOG() << " now " << checklik.LOG();
			errMsg << "Codon61SequenceGibbsSampler::go() loglik diff of " << checklik.LOG()-jlik[i].LOG();
			warning(errMsg.str().c_str());
		}
		rv.change_value(pos,i,Variable::_ACCEPT);
	}
	_mcmc->record_proposal();
	_mcmc->set_alpha(mydouble(1.0));
	_mcmc->set_accept(true);
}
	
} // namespace genomegaMap

