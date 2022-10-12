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
   Takes as input the filename, the number of chains, the length of the chains, integer 0/1 statement (depending on whether the chains are closed, i.e. rings, or open, i.e. linear, 0: ring, 1: linear) and the box dimension (optional)
   If the box dimension is not specified, or if it is equal to 0, then the system is not periodic and the coordinates are unwrapped. 
   If a non-zero box-dimension is specified, the coordinates are unwrapped, according to the PBC.
   Returns the writhe of each chain
   Last line returns the mean absolute writhe over all chains
------------------------------------------------------------------------- */

using namespace std;

int main(int argc, char* argv[])
{
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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
        int chunk = num_chains / size;
        int offset = chunk * rank;
        bool last = false;
        if (rank + 1 == size)
                last = true;
        double result[chunk];
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
        string outfile_name = to_string(chain_length) + "_wr_mpi_out_" + to_string(rank) + ".txt";
        ofstream outfile;
        if (!fs::exists("./output/wr_mpi"))
        {
                cout << "Creating wr_mpi directory..." << endl;
                fs::create_directory("./output/wr_mpi");
        }
        outfile.open("./output/wr_mpi/" + outfile_name);
	//ofstream outfile;
	int count = 0;
	double sum = 0;
	//outfile.open("./output/wr_mpi"+ outfile_name);
        if (ringlinear==0)
        {
           is_closed = true;
        }
        else
        { 
           is_closed = false;
        }
	//for (int i = 0; i < num_chains; i++)
        for (int i = rank * chunk; i < (rank + 1) * chunk; i++)
	{
		double** chain1 = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			chain1[j] = new double[3];
			chain1[j][0] = coords[j + (i * chain_length)][0];
			chain1[j][1] = coords[j + (i * chain_length)][1];
			chain1[j][2] = coords[j + (i * chain_length)][2];
		}

		double res = wr(chain1, chain_length, is_closed);
                //cout << "the result is "<< res <<"\n";
		outfile << res << "\n";;
		sum += abs(res);
		count++;
                result[i - (rank * chunk)] = res;
		delete_array(chain1, chain_length);
	}
        //for (int i=1;i<chunk+1;i++)
        //{
        //   outfile<<result[i-1]<<"\n";
        //}

	delete_array(coords, num_chains);
	outfile.close();
	cout << "Absolute avg wr: " << sum / count << "\n";
        MPI_Finalize();
	return 0;
}
