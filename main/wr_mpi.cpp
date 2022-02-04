#include "../include/funcs.h"
#include "mpi.h"

using namespace std;

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int num_chains = stoi(argv[3]);
	int chain_length = stoi(argv[2]);
	int chunk = num_chains / size;
	int num;
	double** coords = read_coords(argv[1], &num);
	double result[chunk];
	create_ouput_dir();
	string file_name = to_string(chain_length) + "_wr_mpi_out_" + to_string(rank) + ".txt";
	ofstream outfile;
	outfile.open("./output/" + file_name);

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

		double res = wr(chain1, chain_length, false);
		outfile << "writhe of chain " << i << ": " << res << "\n";
		result[i - (rank * chunk)] = res;
		delete_array(chain1, chain_length);
	}

	delete_array(coords, num_chains);
	outfile.close();

	MPI_Finalize();
	return 0;
}
