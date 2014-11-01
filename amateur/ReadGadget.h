/***************************************************************************
 *   Copyright (C) 2010 by Pascal Stephan Philipp Steger                   *
 *   psteger@phys.ethz.ch                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, ma  02111-1307, USA.             *
 ***************************************************************************/

/*!
 * \file ReadGadget.h
 * \brief All publicly known functions from ReadGadget.c and template instantiations.
 */

#ifndef READGADGET_H
#define READGADGET_H

#include <string>

extern void swap_Nbyte(char* data, int n, int m);
size_t my_fread(void* ptr, size_t size, size_t nmemb, FILE* stream);
extern int find_block(FILE* fd, char* label);

template< typename T >
int read_gadget_1( T *data, std::string label, FILE* fd );
//extern int read_gadget_1( float* data, char* label, FILE* fd);
//extern int read_gadget_1( int* data, char* label,FILE* fd);
//extern int read_gadget_1( unsigned* data, char* label, FILE* fd);

template< typename T >
int read_gadget_3( T* data, std::string label, FILE* fd );
//extern int read_gadget_3( float* data, char* label, FILE* fd);

extern int read_gadget_head( int* npart, double* massarr, double* time, double* redshift, FILE* fd );
extern void get_particle_number ( std::string filename, unsigned * id_min, unsigned * id_max, double * masses, unsigned* npart, unsigned& ntot , double& redshift, double& boxsize );

#endif
