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

/**
 * \file Exceptions.h
 * \brief Contains the possible Exceptions that can be thrown.
 */

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>

/**
 * \brief Exception: No halo distribution computed or found
 * This exception is called if operations are attempted on the halo distribution,
 * but no such halo distribution has been calculated (through FindHalos()) or
 * loaded (not yet available).
 */
class ErrNoHalosFound : public std::runtime_error {
 public:
  ErrNoHalosFound() : std::runtime_error("Error : No halo distribution computed or found.") { }
};


/** 
 * \brief Exception: Maximum number of halo progenitors exceeded
 * This exception is called if a progenitor halo is added to a given halo but that
 * halo already possesses the maximum allowed number of progenitor haloes (standard 25).
 * This happens for memory reasons and means you have to increase the number of 
 * allowed progenitors.
 */
class MaxProgenitorsExceed : public std::runtime_error {
 public:
  MaxProgenitorsExceed() : std::runtime_error("Error : Maximum number of progenitor halos exceeded.") {}
};

/**
 * \brief Exception: Wrong dimension (n!=3) in cross product computatoin
 * This exception is called whenever a vector has a wrong dimension when it should be used in a cross product
 * computation.
 */
class CrossProductWrongDimension : public std::runtime_error {
 public:
  CrossProductWrongDimension() : std::runtime_error("Error : Wrong dimension (N!=3) in cross product computation.") {}
};

#endif
