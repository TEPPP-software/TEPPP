#include "../include/funcs.h"
#include "mpi.h"

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
	double** coords = read_coords(argv[1], &num);
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
