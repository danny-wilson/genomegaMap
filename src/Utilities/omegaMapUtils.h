/*  Copyright 2018 Daniel Wilson.
 *
 *  Part of the omegaMap library.
 *
 *  The omegaMap library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The omegaMap library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the omegaMap library. If not, see <http://www.gnu.org/licenses/>.
 */
/*
 *  omegaMapUtils.h
 *  gcat
 *
 *  Created by Daniel Wilson on 10/15/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _OMEGAMAPUTILS_H_
#define _OMEGAMAPUTILS_H_

namespace omegaMapUtils {
	// Variables
	extern const double SQRT2PI;
	extern const char* oneLetterCodes;

	// Functions
	double LambertW(const double z);
	double LambertW1(const double z);
	double hypergeometric1F1(const double a, const double b, const double c);
};

#endif // _OMEGAMAPUTILS_H_
