/****************************************************************************
 * Copyright (c) 2007 Einir Valdimarsson and Chrysanthe Preza
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 ****************************************************************************/
 
#ifndef _ESTIMATE_IO_H
#define _ESTIMATE_IO_H
 
#include <blitz/array.h>

namespace cosm {
 
template<typename T, int N> 
class EstimateIO {
public:
    EstimateIO(){};
    virtual ~EstimateIO() {};
    virtual void ReadData(blitz::Array<T,N>& img, const std::string& name ) {};
    virtual void WriteData(blitz::Array<T,N>& img, const std::string& name ) {};
    virtual void ReadData(blitz::Array<std::complex<T>,N>& img, const std::string& name ) {};
    virtual void WriteData(blitz::Array<std::complex<T>,N>& img, const std::string& name ) {};
};

};
#endif // _ESTIMATE_IO
