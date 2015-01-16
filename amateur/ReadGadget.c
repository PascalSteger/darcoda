/***************************************************************************
*   readgadget.c modified by Pascal Stephan Philipp Steger                 *
*   psteger@phys.ethz.ch                                                   *
****************************************************************************/

/*!
 * \file ReadGadget.c
 * \brief sample input file from Gadget2, modified for use with AMATEUR
 */

#include <iostream>
#include <string>
#include <string.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>

// for low level routines from read_gadget.c
int blksize, swap_rg = 0;
#define SKIP  {my_fread(&blksize,sizeof(int),1,fd); swap_Nbyte((char*)&blksize,1,sizeof(int));}


/* \brief Basic routine to read data from a file */
/* Reads a slab from a given filestream
 * \param ptr: pointer for fread
 * \param size: size of slab
 * \param nmemb:
 * \param stream: filestream
 * \return nread: how many kilobytes were read?
 */
size_t my_fread ( void * ptr, size_t size, size_t nmemb, FILE * stream ) {
	size_t nread;

	if ( ( nread = fread ( ptr, size, nmemb, stream ) ) != nmemb ) {
		std::cerr << "I/O error (fread) !" << std::endl;
		exit ( 3 );
	}
	return nread;
}

/**
 * swap endian
 * @param data
 * @param n how many entries do we have in array?
 * @param m how many bytes are there for a bigger number type
 */
void swap_Nbyte( char *data, unsigned n, unsigned m ) {
	char old_data[ 16 ];

	if ( swap_rg > 0 ) {
		for ( unsigned i = 0; i < n; ++i ) {
			memcpy( &old_data[ 0 ], &data[ i * m ], m );
			for ( unsigned j = 0; j < m; ++j ) {
				data[ i * m + j ] = old_data[ m - j - 1 ];
			}
		}
	}
}

/**
 * find a block in snapshot file
 * @param fd
 * @param label
 * @return
 */
int find_block ( FILE *fd, std::string label ) {
	int blocksize = 0;
	char label2[5] = {"    "};
	label.copy(label2,4);
	char blocklabel[ 5 ] = {"    "};

	rewind ( fd );

	while ( !feof ( fd ) && blocksize == 0 ) {
		SKIP;
		if ( blksize == 134217728 ) {
			std::cerr << "Enable ENDIAN swapping!" << std::endl;
			swap_rg = 1 - swap_rg;
			swap_Nbyte ( ( char* ) & blksize, 1, sizeof( int ) );
		}
		if ( blksize != 8 ) {
			std::cerr << "Incorrect format (blksize = " << blksize << "d)!" << std::endl;
			exit ( 1 );
		} else {
			my_fread ( blocklabel, 4 * sizeof ( char ), 1, fd );
			my_fread ( &blocksize, sizeof ( int ), 1, fd );
			swap_Nbyte ( ( char* ) & blocksize, 1, 4 );
			// 			std::cout << "Found Block " << blocklabel
			// 			<< " with " << blocksize << " bytes." << std::endl;
			SKIP;
			if ( strcmp ( label2, blocklabel ) != 0 ) {
				fseek ( fd, blocksize, 1 );
				blocksize = 0;
			}
		}
	}
	return ( blocksize -8 );
}

/**
 * read gadget header info
 * \param npart
 * \param Massarr
 * \param time
 * \param redshift
 * \param FlagSfr
 * \param FlagFeedback
 * \param Nall
 * \param FlagCooling
 * \param NumFiles
 * \param BoxSize
 * \param fd
 * \return size of header block
 */
int read_gadget_head ( unsigned* npart, double* Massarr, double* time, double* redshift, int* FlagSfr, int* FlagFeedback, unsigned* Nall, int* FlagCooling, int* NumFiles, double* BoxSize, FILE* fd ) {
	int blocksize, dummysize;

	blocksize = find_block ( fd, "HEAD" );
	if ( blocksize <= 0 ) {
		std::cerr << "Block " << "HEAD" << " not fond!" << std::endl;
		exit ( 5 );
	}

	dummysize = blocksize - 6 * sizeof ( int ) - 8 * sizeof ( double );
	SKIP;
	my_fread ( npart, 6 * sizeof ( unsigned ), 1, fd );
	swap_Nbyte ( ( char* ) npart, 6, sizeof( unsigned ) );
	my_fread ( Massarr, 6 * sizeof ( double ), 1, fd );
	swap_Nbyte ( ( char* ) Massarr, 6, sizeof( double ) );
	my_fread ( time, sizeof ( double ), 1, fd );
	swap_Nbyte ( ( char* ) time, 1, sizeof( double ) );
	my_fread ( redshift, sizeof ( double ), 1, fd );
	swap_Nbyte ( ( char* ) redshift, 1, sizeof( double ) );
	my_fread ( FlagSfr, sizeof( int ), 1, fd );
	swap_Nbyte ( ( char* ) FlagSfr, 1, sizeof( int ) );
	my_fread ( FlagFeedback, sizeof( int ), 1, fd );
	swap_Nbyte ( ( char* ) FlagFeedback, 1, sizeof( int ) );
	my_fread ( Nall, 6 * sizeof ( int ), 1, fd );
	swap_Nbyte ( ( char* ) Nall, 6, sizeof( int ) );
	my_fread ( FlagCooling, sizeof( int ), 1, fd );
	swap_Nbyte ( ( char* ) FlagCooling, 1, sizeof( int ) );
	my_fread ( NumFiles, sizeof( int ), 1, fd );
	swap_Nbyte ( ( char* ) NumFiles, 1, sizeof( int ) );
	my_fread ( BoxSize, sizeof ( double ), 1, fd );
	swap_Nbyte ( ( char* ) BoxSize, 1, sizeof( double ) );
	fseek ( fd, dummysize, 1 );
	SKIP;

	return blocksize;
}


/**
 * read a 1D array
 * @param data
 * @param label
 * @param fd
 * @return
 */
template <typename T>
int read_gadget_1 ( T *data, std::string label, FILE *fd ) {
	int blocksize;
	char label2[5] = {"    "};
	label.copy(label2,4);

	blocksize = find_block ( fd, label2 );
	if ( blocksize <= 0 ) {
		std::cerr << "Block " << label2 << " not found!" << std::endl;
		exit ( 5 );
	} else {
		std::cout << "Reading " << blocksize << " bytes of data from "
		<< label << "..." << std::endl;
		SKIP;
		my_fread ( data, blocksize, 1, fd );
		swap_Nbyte ( ( char* ) data, blocksize / sizeof ( T ), sizeof( T ) );
		SKIP;
	}
	return blocksize / sizeof ( T );
}


/**
 * read a 3D array
 * @param data
 * @param label
 * @param fd
 * @return
 */
template <typename T>
int read_gadget_3 ( T* data, std::string label, FILE* fd ) {
	int blocksize;
	char label2[] = {"    "};
	label.copy(label2, 4);

	blocksize = find_block ( fd, label );
	if ( blocksize <= 0 ) {
		std::cerr << "Block " << label << " not found!" << std::endl;
		exit ( 5 );
	} else {
		std::cout << "Reading " << blocksize << " bytes of data from " << label << "..." << std::endl;
		SKIP;
		my_fread ( data, blocksize, 1, fd );
		swap_Nbyte ( ( char* ) data, blocksize / sizeof ( T ), sizeof( T ) );
		SKIP;
	}
	return ( blocksize / sizeof ( T ) / 3 );
}

/* \brief Get particlenumber from a snapshot file */
/* Get the number of particles stored inside a snapshot file.
 * \param filename
 * \param id_min
 * \param id_max
 * \param masses
 * \param npart
 * \param ntot
 * \param redshift
 * \param boxsize
 */
void get_particle_number ( std::string filename, unsigned * id_min, unsigned * id_max, double * masses, unsigned * npart, unsigned& ntot , double& redshift, double& boxsize ) {
	FILE * fd = 0;
	int FlagSfr, FlagFeedback, FlagCooling, NumFiles;
	unsigned Nall[ 6 ];
	double time;
	unsigned * id_block;

	if ( ! ( fd = fopen ( filename.c_str(), "r" ) ) ) {
		std::cerr << "Cannot open file " << filename << "!" << std::endl;
		exit ( 2 );
	}

	// read header to get global properties
	n = read_gadget_head ( npart, masses, &time, &redshift, &FlagSfr, &FlagFeedback, Nall, &FlagCooling, &NumFiles, &boxsize, ( FILE * ) fd );

	// show details
	ntot = 0;
	for ( int i = 0; i < 6; i++ ) {
		std::cout << "particle type ";
		std::cout.width( 5 );
		std::cout << i << ", count = ";
		std::cout.width( 10 );
		std::cout << npart[ i ] << ", fixed mass = ";
		std::cout.width( 10 );
		std::cout << masses[ i ] << std::endl;
		ntot += npart[ i ];
	}
	std::cout << "Time of snapshot = " << time << ", z = " << redshift << ", ntot = " << ntot << ", boxsize = " << boxsize << std::endl << std::endl;

	id_block = ( unsigned * ) malloc ( ntot * sizeof ( int ) );

	// read data blocks
	n = read_gadget_1 ( ( unsigned * ) id_block, "ID  ", fd );
	//for ( unsigned i = 0; i < ntot; ++i )
	//	id_block[ i ] = unsigned( id_block[ i ] );
	// search minima and maxima in ID-list
	// use counter to follow the list of particle numbers successively
	unsigned cntr = 0;
	for ( unsigned j = 0; j < 4; ++j ) {
		id_min[ j ] = id_block[ cntr ];
		id_max[ j ] = id_block[ cntr ];
		for ( unsigned i = cntr; i < cntr + npart[ j ] - 1; ++i ) {
			if ( id_min[ j ] > id_block[ i ] )
				id_min[ j ] = id_block[ i ];
			if ( id_max[ j ] < id_block[ i ] )
				id_max[ j ] = id_block[ i ];
		}
		cntr += npart[ j ];
	}
	id_min[ 4 ] = id_max[ 3 ] + 1;
	id_max[ 4 ] = id_min[ 4 ] + npart[ 4 ];
	id_min[ 5 ] = 0;
	id_max[ 5 ] = 0;
	for ( unsigned j = 0; j < 6; ++j ) {
		// give out ID range
		std::cout << "ID(" << j << ") in range: ";
		std::cout.width( 10 );
		std::cout << id_min[ j ] << ",";
		std::cout.width ( 10 );
		std::cout << id_max[ j ] << std::endl;
	}
	std::cout << std::endl;

	fclose ( fd );
	free ( id_block );
}

template int read_gadget_1<float>( float *data, std::string label, FILE *fd );
template int read_gadget_1<unsigned>( unsigned *data,std::string label, FILE *fd );
template int read_gadget_1<int>( int *data,std::string label, FILE *fd );
template int read_gadget_3<float>( float *data, std::string label, FILE *fd );
