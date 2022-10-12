/* -*- -*- ----------------------------------------------------------
   TEPPP: Topological Entanglement in Polymers, Proteins and Periodic structures
   https://github.com/TEPPP-software/TEPPP.git
   Eleni Panagiotou, epanagio@asu.edu

   Copyright (2021) Eleni Panagiotou This software is distributed under
   the BSD 3-Clause License.

   See the README file in the top-level TEPPP directory.
   Contributors: Tom Herschberg, Kyle Pifer and Eleni Panagiotou
------------------------------------------------------------------------- */

#include "../include/funcs.h"
using namespace std;


/* -*- -*- ----------------------------------------------------------
   Takes as input the filename, the number of chains, the length of the chains, integer 0/1 statement (depending on whether the chains are closed, i.e. rings, or open, i.e. linear, 0: ring, 1: linear) and the box dimension (optional)
   If the box dimension is not specified, or if it is equal to 0, then the system is not periodic and the coordinates are unwrapped. 
   If a non-zero box-dimension is specified, the coordinates are unwrapped, according to the PBC.
   Returns the Gauss linking integral between each pair of chains
   Last line returns the mean absolute linking number over all pairs
------------------------------------------------------------------------- */


int main(int argc, char* argv[])
{
	int num_chains = stoi(argv[3]);
	int chain_length = stoi(argv[2]);
	double box_dim;
        bool is_closed;
        int ringlinear = stoi(argv[4]);
	if (argc >= 6)
		box_dim = stod(argv[5]);
	else
		box_dim = 0;
	int num;
	double** coords;
	if (box_dim == 0)
	{
		coords = read_coords(argv[1], &num);
	}
	else
	{
		coords = read_coords(argv[1], &num, chain_length, box_dim);
	}

	ofstream outfile;
	int count = 0;
	double sum = 0;
	create_output_dir();
	outfile.open("./output/lk_out.txt");
        if (ringlinear==0)
        {
           is_closed = true;
        }
        else
        {
           is_closed = false;
        }

	for (int i = 0; i < num_chains; i++)
	{
		double** chain1 = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			chain1[j] = new double[3];
			chain1[j][0] = coords[j + (i * chain_length)][0];
			chain1[j][1] = coords[j + (i * chain_length)][1];
			chain1[j][2] = coords[j + (i * chain_length)][2];
		}

		/**
		 * @brief This for loop is used to calculate the coordinates
		 * for both chain 1 and chain 2 for the lk calculation.
		 *
		 */
		for (int j = i + 1; j < num_chains; j++)
		{
			double** chain2 = new double*[chain_length];
			for (int k = 0; k < chain_length; k++)
			{
				chain2[k] = new double[3];
				chain2[k][0] = coords[k + (j * chain_length)][0];
				chain2[k][1] = coords[k + (j * chain_length)][1];
				chain2[k][2] = coords[k + (j * chain_length)][2];
			}

			double res = lk(chain1, chain2, chain_length, chain_length, is_closed);
			outfile << res << "\n";
			sum += abs(res);
			count++;
			delete_array(chain2, chain_length);
		}

		delete_array(chain1, chain_length);
	}

	delete_array(coords, num_chains);
	outfile.close();
	cout << "Absolute avg lk: " << sum / count << "\n";

	return 0;
}
