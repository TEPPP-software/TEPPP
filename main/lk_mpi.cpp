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
	double result[chunk][num_chains];
	create_output_dir();
	string file_name = to_string(chain_length) + "_lk_mpi_out_" + to_string(rank) + ".txt";
	ofstream outfile;
	outfile.open("./output/" + file_name);

	for (int i = rank * chunk; i < (rank + 1) * chunk; i++)
	{
		result[i - (rank * chunk)][i] = 0;
		double** chain1 = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			chain1[j] = new double[3];
			chain1[j][0] = coords[j + (i * chain_length)][0];
			chain1[j][1] = coords[j + (i * chain_length)][1];
			chain1[j][2] = coords[j + (i * chain_length)][2];
		}

		for (int j = 0; j < num_chains; j++)
		{
			if (i == j)
				continue;
			double** chain2 = new double*[chain_length];
			for (int k = 0; k < chain_length; k++)
			{
				chain2[k] = new double[3];
				chain2[k][0] = coords[k + (j * chain_length)][0];
				chain2[k][1] = coords[k + (j * chain_length)][1];
				chain2[k][2] = coords[k + (j * chain_length)][2];
			}

			double res = lk(chain1, chain2, chain_length, chain_length, false);
			outfile << "linking number between chains " << i << " and " << j << ": " << res << "\n";
			result[i - (rank * chunk)][j] = res;
			delete_array(chain2, chain_length);
		}

		delete_array(chain1, chain_length);
	}

	delete_array(coords, num_chains);
	outfile.close();

	MPI_Finalize();
	return 0;
}
