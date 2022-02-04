/*
	read_dcd : c++ class + main file example for reading a CHARMM dcd file
	Copyright (C) 2013  Florent Hedin

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "../include/dcd_r.hpp"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
using namespace std;

DCD_R::DCD_R(const char filename[])
{
	dcdf.exceptions(std::ifstream::failbit);
	try
	{
		dcdf.open(filename, ios::in | ios::binary);
	}
	catch (std::ifstream::failure e)
	{
		cerr << "Exception opening/reading file '" << filename << "' : " << std::endl;
		cerr << "Please chech the path of the file and if it exists." << endl;
	}

	dcd_first_read = true;
}

void DCD_R::alloc()
{
	X = new float[NATOM];
	Y = new float[NATOM];
	Z = new float[NATOM];
	pbc[0] = pbc[1] = pbc[2] = pbc[3] = pbc[4] = pbc[5] = 0.0;
}

void DCD_R::read_header()
{
	unsigned int fortcheck1, fortcheck2;

	// This is the trick for reading binary data from fortran file : see the method
	// DCD::checkFortranIOerror for more details. we are reading data corresponding to a "write(...)
	// HDR,ICNTRL" fortran statement
	dcdf.read((char*) &fortcheck1, sizeof(unsigned int)); // consistency check 1
	dcdf.read((char*) HDR,
			  sizeof(char) * 4); // first data block written by fortran  : a character array of length 4.
	dcdf.read(
		(char*) ICNTRL,
		sizeof(int) * 20); // second data block written by fortran : an integer(4) array of length 20.
	dcdf.read((char*) &fortcheck2, sizeof(unsigned int)); // consistency check 2
	checkFortranIOerror(
		__FILE__,
		__LINE__,
		fortcheck1,
		fortcheck2); // if the 2 unsigned ints have a different value there was an error

	/* See dcd.hpp for details on ICNTRL */
	HDR[4] = '\0';
	NFILE = ICNTRL[0];
	NPRIV = ICNTRL[1];
	NSAVC = ICNTRL[2];
	NSTEP = ICNTRL[3];
	NDEGF = ICNTRL[7];
	FROZAT = ICNTRL[8];
	DELTA4 = ICNTRL[9];
	QCRYS = ICNTRL[10];
	CHARMV = ICNTRL[19];

	/* Several "lines" of title of length 80 are written to the dcd file by CHARMM */
	dcdf.read((char*) &fortcheck1, sizeof(unsigned int));
	dcdf.read((char*) &NTITLE, sizeof(int));
	if (NTITLE == 0)
	{
		TITLE = new char[80 + 1];
		TITLE[0] = '\0';
	}
	else
	{
		TITLE = new char[NTITLE * 80 + 1];
		for (int it = 0; it < NTITLE; it++)
		{
			dcdf.read((char*) &TITLE[it * 80], sizeof(char) * 80);
		}
		TITLE[NTITLE * 80] = '\0';
	}
	dcdf.read((char*) &fortcheck2, sizeof(unsigned int));
	checkFortranIOerror(__FILE__, __LINE__, fortcheck1, fortcheck2);

	// reading number of atoms
	dcdf.read((char*) &fortcheck1, sizeof(unsigned int));
	dcdf.read((char*) &NATOM, sizeof(int));
	dcdf.read((char*) &fortcheck2, sizeof(unsigned int));
	checkFortranIOerror(__FILE__, __LINE__, fortcheck1, fortcheck2);

	/* If some atoms of the MD or MC simulation are frozen (i.e. never moving ) it is useless to store
	 * their coordinates more than once. In that case a list of Free atoms (moving ones) is written at
	 * the end of the header part of the dcd. See DCD_R::read_oneFrame() for more details.
	 */
	LNFREAT = NATOM - FROZAT;
	if (LNFREAT != NATOM)
	{
		FREEAT = new int[LNFREAT];
		dcdf.read((char*) &fortcheck1, sizeof(unsigned int));
		dcdf.read((char*) FREEAT, sizeof(int) * LNFREAT);
		dcdf.read((char*) &fortcheck2, sizeof(unsigned int));
		checkFortranIOerror(__FILE__, __LINE__, fortcheck1, fortcheck2);
	}

	// allocate memory for storing coordinates (only one frame of the dcd is stored, so several
	// (NFILE) calls to DCD_R::read_oneFrame() are necessary for reading the whole file).
	alloc();
}

void DCD_R::read_oneFrame()
{
	unsigned int fortcheck1, fortcheck2;

	int siz = (dcd_first_read) ? NATOM : LNFREAT;

	float* tmpX = new float[siz];
	float* tmpY = new float[siz];
	float* tmpZ = new float[siz];

	if (QCRYS)
	{
		dcdf.read((char*) &fortcheck1, sizeof(unsigned int));
		dcdf.read((char*) pbc, sizeof(double) * 6);
		dcdf.read((char*) &fortcheck2, sizeof(unsigned int));
		checkFortranIOerror(__FILE__, __LINE__, fortcheck1, fortcheck2);
	}

	// X
	dcdf.read((char*) &fortcheck1, sizeof(unsigned int));
	dcdf.read((char*) tmpX, sizeof(float) * siz);
	dcdf.read((char*) &fortcheck2, sizeof(unsigned int));
	checkFortranIOerror(__FILE__, __LINE__, fortcheck1, fortcheck2);

	// Y
	dcdf.read((char*) &fortcheck1, sizeof(unsigned int));
	dcdf.read((char*) tmpY, sizeof(float) * siz);
	dcdf.read((char*) &fortcheck2, sizeof(unsigned int));
	checkFortranIOerror(__FILE__, __LINE__, fortcheck1, fortcheck2);

	// Z
	dcdf.read((char*) &fortcheck1, sizeof(unsigned int));
	dcdf.read((char*) tmpZ, sizeof(float) * siz);
	dcdf.read((char*) &fortcheck2, sizeof(unsigned int));
	checkFortranIOerror(__FILE__, __LINE__, fortcheck1, fortcheck2);

	if (dcd_first_read)
	{
		memcpy(X, tmpX, NATOM * sizeof(float));
		memcpy(Y, tmpY, NATOM * sizeof(float));
		memcpy(Z, tmpZ, NATOM * sizeof(float));
	}
	else
	{
		if (LNFREAT != NATOM)
		{
			for (int it = 0; it < siz; it++)
			{
				X[FREEAT[it] - 1] = tmpX[it];
				Y[FREEAT[it] - 1] = tmpY[it];
				Z[FREEAT[it] - 1] = tmpZ[it];
			}
		}
		else
		{
			memcpy(X, tmpX, NATOM * sizeof(float));
			memcpy(Y, tmpY, NATOM * sizeof(float));
			memcpy(Z, tmpZ, NATOM * sizeof(float));
		}
	}

	if (dcd_first_read)
		dcd_first_read = false;

	delete[] tmpX;
	delete[] tmpY;
	delete[] tmpZ;
}

void DCD_R::printHeader() const
{
	int i;

	cout << "HDR :\t" << HDR << endl;

	cout << "ICNTRL :\t";
	for (i = 0; i < 20; i++)
		cout << ICNTRL[i] << "\t";
	cout << endl;

	cout << "NTITLE :\t" << NTITLE << endl;
	cout << "TITLE :\t" << TITLE << endl;

	cout << "NATOM :\t" << NATOM << endl;
	cout << "LNFREAT :\t" << LNFREAT << endl;
}

DCD_R::~DCD_R()
{
	dcdf.close();

	delete[] TITLE;

	if (LNFREAT != NATOM)
		delete[] FREEAT;

	delete[] X;
	delete[] Y;
	delete[] Z;
}
