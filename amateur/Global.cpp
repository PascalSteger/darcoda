/**************************************************************************
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
 * \file Global.cpp
 * \brief Implementation of all system-wide properties that all parts of the program should know.
 */

#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <limits>
#include <math.h>
#include "Global.h"


/*
 * checks whether a file does exist
 * @param filename
 * @return bool true or false
 */
bool DoesFileExist( char *filename ) {
  bool flag = false;
  std::fstream fin( filename, std::ios::in );
  if ( fin.is_open() )
    flag = true;
  fin.close();
  return flag;
}

/**
 *\brief print characteristics of basic types
 * Show Minimum, Maximum, Sign, Bits, Inf for bool, char, int, long, unsigned, unsigned long, float, double
 *\return: void
 */
void
printMinMaxValues(){

  std::cout << "Minimum value for bool: " << std::numeric_limits<bool>::min() << std::endl;
  std::cout << "Maximum value for bool: " << std::numeric_limits<bool>::max() << std::endl;
  std::cout << "bool is signed: " << std::numeric_limits<bool>::is_signed << std::endl;
  std::cout << "Non-sign bits in bool: " << std::numeric_limits<bool>::digits << std::endl;
  std::cout << "bool has infinity: " << std::numeric_limits<bool>::has_infinity << std::endl;

  std::cout << "Minimum value for char: " << std::numeric_limits<char>::min() << std::endl;
  std::cout << "Maximum value for char: " << std::numeric_limits<char>::max() << std::endl;
  std::cout << "char is signed: " << std::numeric_limits<char>::is_signed << std::endl;
  std::cout << "Non-sign bits in char: " << std::numeric_limits<char>::digits << std::endl;
  std::cout << "char has infinity: " << std::numeric_limits<char>::has_infinity << std::endl;

  std::cout << "Minimum value for int: " << std::numeric_limits<int>::min() << std::endl;
  std::cout << "Maximum value for int: " << std::numeric_limits<int>::max() << std::endl;
  std::cout << "int is signed: " << std::numeric_limits<int>::is_signed << std::endl;
  std::cout << "Non-sign bits in int: " << std::numeric_limits<int>::digits << std::endl;
  std::cout << "int has infinity: " << std::numeric_limits<int>::has_infinity << std::endl;

  std::cout << "Minimum value for long: " << std::numeric_limits<long>::min() << std::endl;
  std::cout << "Maximum value for long: " << std::numeric_limits<long>::max() << std::endl;
  std::cout << "long is signed: " << std::numeric_limits<long>::is_signed << std::endl;
  std::cout << "Non-sign bits in long: " << std::numeric_limits<long>::digits << std::endl;
  std::cout << "long has infinity: " << std::numeric_limits<long>::has_infinity << std::endl;

  std::cout << "Minimum value for unsigned: " << std::numeric_limits<unsigned>::min() << std::endl;
  std::cout << "Maximum value for unsigned: " << std::numeric_limits<unsigned>::max() << std::endl;
  std::cout << "unsigned is signed: " << std::numeric_limits<unsigned>::is_signed << std::endl;
  std::cout << "Non-sign bits in unsigned: " << std::numeric_limits<unsigned>::digits << std::endl;
  std::cout << "unsigned has infinity: " << std::numeric_limits<unsigned>::has_infinity << std::endl;

  std::cout << "Minimum value for unsigned long: " << std::numeric_limits<unsigned long>::min() << std::endl;
  std::cout << "Maximum value for unsigned long: " << std::numeric_limits<unsigned long>::max() << std::endl;
  std::cout << "unsigned long is signed: " << std::numeric_limits<unsigned long>::is_signed << std::endl;
  std::cout << "Non-sign bits in unsigned long: " << std::numeric_limits<unsigned long >::digits << std::endl;
  std::cout << "unsigned long has infinity: " << std::numeric_limits<unsigned long>::has_infinity << std::endl;

  std::cout << "Minimum value for float: " << std::numeric_limits<float>::min() << std::endl;
  std::cout << "Maximum value for float: " << std::numeric_limits<float>::max() << std::endl;
  std::cout << "float is signed: " << std::numeric_limits<float>::is_signed << std::endl;
  std::cout << "Non-sign bits in float: " << std::numeric_limits<float>::digits << std::endl;
  std::cout << "float has infinity: " << std::numeric_limits<float>::has_infinity << std::endl;

  std::cout << "Minimum value for double: " << std::numeric_limits<double>::min() << std::endl;
  std::cout << "Maximum value for double: " << std::numeric_limits<double>::max() << std::endl;
  std::cout << "double is signed: " << std::numeric_limits<double>::is_signed << std::endl;
  std::cout << "Non-sign bits in double: " << std::numeric_limits<double>::digits << std::endl;
  std::cout << "double has infinity: " << std::numeric_limits<double>::has_infinity << std::endl;
}

/*
 * write a given string, without newline
 * @param text: this string
 */
void write(std::string text){
  std::cout << text;
  std::cout << std::flush;
}

/* \brief output a comment on a new line */
/* Write a string onto a line on std::out and add a newline
 * @param text: comment to be printed
 */
void writeln(std::string text){
  std::cout << text << std::endl;
  std::cout << std::flush;
}

/* 
 * Write "done" without newline
 */
void writedone(void){
  std::cout << "done!";
  std::cout << std::flush;
}

/*
 * Write "done" with newline
 */
void writelndone(void){
  std::cout << "done!" << std::endl;
  std::cout << std::flush;
}

/*
 * Write "done!" if DebugLevel activated
 */
void donedebug(bool DebugLevelActivated){
  if(DebugLevelActivated)
    writelndone();
}

/*
 * Write a string onto a line on std::cout,
 * given that the first argument is true
 * @param DebugLevelActivated: boolean value determining output funcionality
 * @param text: std::string of desired output
 */
void writedebug(bool DebugLevelActivated, std::string text){
  if(DebugLevelActivated)
    write(text);
}

/* \brief output a debug comment on std out */
/* Write a string onto a line on std::out and add a newline,
 * given that the first argument is true.
 * @param DebugLevelActivated: boolean value determining output functionality
 * @param text: std::string of desired output
 */
void writelndebug(bool DebugLevelActivated, std::string text){
  if(DebugLevelActivated)
    writeln(text);
}

/* \brief run a given command */
/* Run a command with system()
 * @param command: command that should be executed
 */
void RunCommand( std::string command ){
  system(command.c_str());
}

/* @brief Count the lines of a textfile */
/**
 * count the lines of a given file
 * @param filename: path to the file
 * @return unsigned number of lines
 */
unsigned GetLineNumber(std::string filename){
  std::ifstream infile(filename.c_str());
  if(!infile)
    return 0;
  unsigned i = 0;
  std::string sd;
  while(getline(infile,sd) && sd.size()>0)
    ++i;

  infile.close();
  return i;
}

/* @brief Reverse a string */
/**
 * Reverse the order of characters in a string
 * @param s: array of chars (<=> string)
 * @return void
 */
void reverse(char s[])
{
  int c, i, j;

  for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
    c = s[i];
    s[i] = s[j];
    s[j] = c;
  }
}

/* \brief Conversion int -> char[] */
/**
 * Converts an integer into a char[ ].
 * @param n: integer that should be printed as a character array
 * @param s: char array that should hold the number
 * @return void
 */
void itoa(int n, char s[])
{
  int i, sign;

  if ((sign = n) < 0)  /* record sign */
    n = -n;          /* make n positive */
  i = 0;
  do {       /* generate digits in reverse order */
    s[i++] = n % 10 + '0';   /* get next digit */
  } while ((n /= 10) > 0);     /* delete it */
  if (sign < 0)
    s[i++] = '-';
  s[i] = '\0';
  reverse(s);
}

/* \brief Conversion int -> std::string */
/**
 * Convert an integer into a string.
 * @param intValue: integer to be converted
 * @return std::string of that integer value
 */
std::string IntToString(int intValue) {
  char *myBuff;
  std::string strRetVal;
  myBuff = new char[100];
  memset(myBuff,'\0',100);

  itoa(intValue,myBuff);

  strRetVal = myBuff;
  delete[] myBuff;
  return(strRetVal);
}

/* \brief returns folder position */
/**
 * Returns the paths to all important folders.
 * @param where: which folder? allowed: amd, snap, amateur, output, parameters
 * @return std::string full path
 */
std::string Folder(std::string where){
  std::string folder_amd="/data/achtland1/psteger/amd/";
  std::string folder_snaps = folder_amd+"snap/";
  std::string folder_amateur = folder_amd+"amateur/";
  std::string folder_output = folder_amd+"output/";
  std::string folder_parameters = folder_output + "parameters/";
  if(where=="amd")
    return folder_amd;
  else if(where=="snaps")
    return folder_snaps;
  else if(where=="amateur")
    return folder_amateur;
  else if(where=="output")
    return folder_output;
  else if(where=="parameters")
    return folder_parameters;
  else
    {
      std::cout << "Folder " << where << " not found, returning folder of amd!";
      return folder_amd;
    }
}

/* \brief Returns string of an AMATEUR snapshot */
/**
 * Returns the full path to an AMATEUR snapshot
 * @param SnapNo: ID of snapshot
 * @return std::string full_path
 */
std::string Snaps(unsigned SnapNo){
  std::string snaps [18]= {

    "tarkin_5726_10x_snap_302", 
    "tarkin_25174_CSF_10x_snap_008", 
    "tarkin_25174_CSF_10x_snap_009",
    "tarkin_13309_CSF_45x_snap_002",
    "tarkin_13309_CSF_45x_snap_058",
    "tarkin_13309_CSF_45x_snap_118",
    "tarkin_13309_CSF_45x_snap_152",
    "tarkin_13309_CSF_45x_snap_179",
    "tarkin_13309_CSF_45x_snap_200",
    "tarkin_13309_CSF_45x_snap_211",
    "tarkin_13309_CSF_45x_snap_225",
    "tarkin_13309_CSF_45x_snap_257",
    "tarkin_13309_CSF_45x_snap_279",
    "tarkin_13309_CSF_45x_snap_299",
    "tarkin_13309_CSF_45x_snap_302",
    "tarkin_21926_CSF_45x_snap_302",
    "tarkin_25174_CSF_45x_snap_289",
    "tarkin_25174_CSF_45x_snap_302"
  };

  return snaps[SnapNo];
}


real
fbmod( real x, real m ){
  return x - floor( x/m + 0.5 )*m;
}

template< typename real_t >
real_t
SPHK3D( const real_t r, const real_t h ){
  const real_t xi(r/h);
  const real_t norm(8.0/PI/(h*h*h));

  if( xi<0.5 ){
    return norm*(6.0*xi*xi*(xi-1.0)+1.0);
  } else {
    if( xi <= 1.0 ){
      real_t t(1.0-xi);
      return norm*2.0*t*t*t;
    }
    return 0.0;
  }
}
