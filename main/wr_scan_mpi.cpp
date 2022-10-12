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
#include "mpi.h"

/* -*- -*- ----------------------------------------------------------
   Takes as input the filename, the number of chains, the length of the chains, the starting interval, the end interval, the step size and the box dimension (optional)
   If the box dimension is not specified, or if it is equal to 0, then the system is not periodic and the coordinates are unwrapped. 
   If a non-zero box-dimension is specified, the coordinates are unwrapped, according to the PBC.
   Returns the writhe of each interval of a chain

------------------------------------------------------------------------- */


using namespace std;

int main(int argc, char* argv[])
{
	if (argc < 7)
	{
		cout << "Not enough parameters! Exiting...\n";
		return 0;
	}
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int chain_length = stoi(argv[2]);
	int num_chains = stoi(argv[3]);
	int start_chunk = stoi(argv[4]);
	int end_chunk = stoi(argv[5]);
	int step = stoi(argv[6]);
        int chunk = num_chains / size;
	int num;
        double box_dim;
        if (argc >= 7)
                box_dim = stod(argv[6]);
        else
                box_dim = 0;
        double** coords;
        if (box_dim == 0)
        {
                coords = read_coords(argv[1], &num);
        }
        else
        {
                coords = read_coords(argv[1], &num, chain_length, box_dim);
        }
	create_output_dir();
        string file_name = to_string(chain_length) + "_wr_scan_mpi_out_" + to_string(rank) + ".txt";
        if (!fs::exists("./output/wr_scan_mpi"))
        {
                cout << "Creating wr_scan_mpi directory..." << endl;
                fs::create_directory("./output/wr_scan_mpi");
        }
        ofstream outfile;
        outfile.open("./output/wr_scan_mpi/" + file_name);
	//ofstream outfile;
	//outfile.open("./output/wr_scan_out.txt");

	for (int i = rank*chunk; i < (rank+1)*chunk; i++)
	{
		double** chain1 = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			chain1[j] = new double[3];
			chain1[j][0] = coords[j + (i * chain_length)][0];
			chain1[j][1] = coords[j + (i * chain_length)][1];
			chain1[j][2] = coords[j + (i * chain_length)][2];
		}

		for (int a = start_chunk; a <= end_chunk; a += step)
		{
			double** temp_chain1 = new double*[a];
			for (int b = 0; b < a; b++)
			{
				temp_chain1[b] = new double[3];
			}

			for (int b = 0; b < chain_length - a; b++)
			{
				for (int c = 0; c < a; c++)
				{
					temp_chain1[c][0] = chain1[b + c][0];
					temp_chain1[c][1] = chain1[b + c][1];
					temp_chain1[c][2] = chain1[b + c][2];
				}

				double res = wr(temp_chain1, a, false);
				outfile << "writhe of chain " << i << " at atom " << b << " with chunk length " << a << ": " << res
						<< "\n";
			}

			delete_array(temp_chain1, a);
		}

		delete_array(chain1, chain_length);
	}

	delete_array(coords, num_chains);
	outfile.close();
        MPI_Finalize();
	return 0;
}
