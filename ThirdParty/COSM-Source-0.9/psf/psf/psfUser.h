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
                                                                                
#ifndef _PSF_USER_H
#define _PSF_USER_H

#include <iostream>

namespace cosm {

class PsfUser {

public:
    enum Type {
	UNKNOWN = 0x0,
	RADIAL = 0x1,
	COMPLETE = 0x2
    };

public:
    PsfUser() {};
    virtual ~PsfUser() {};
    virtual void update( Type type, int count, int total) 
    { std::cout << (type == RADIAL ? "RADIAL" : "COMPLETE") <<", count: "<< count <<", total: "<< total << std::endl; };

protected:
    // not allowed
    PsfUser(PsfUser&);
    PsfUser& operator=(PsfUser&);
	
};

}

#endif // _PSF_USER_H;
