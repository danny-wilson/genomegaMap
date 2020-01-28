/*  Copyright 2018 Daniel Wilson.
 *
 *  genomegaMapLibrary.h
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
#ifndef _GENOMEGAMAP_LIBRARY_H_
#define _GENOMEGAMAP_LIBRARY_H_
#include <DAG/DAG.h>

// Standard name that gcat looks for when loading the dynamic library
extern "C" gcat::xsd_string load_gcat_library();

#endif // _GENOMEGAMAP_LIBRARY_H_
