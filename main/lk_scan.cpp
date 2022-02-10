#include <filesystem>

#include "../include/funcs.h"
namespace fs = std::filesystem;
using namespace std;

int main(int argc, char* argv[])
{
	if (argc < 7)
	{
		cout << "Not enough parameters! Exiting...\n";
		return -1;
	}

	int chain_length = stoi(argv[2]);
	int num_chains = stoi(argv[3]);
	int start_chunk = stoi(argv[4]);
	int end_chunk = stoi(argv[5]);
	int step = stoi(argv[6]);
	double box_dim;
	if (argc >= 8)
		box_dim = stod(argv[7]);
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
	create_output_dir();
	ofstream outfile;
	outfile.open("./output/lk_scan_out.txt");

	for (int i = 0; i < num_chains; i++)
	{
		/**
		 *  Creating the entire chain1
		 */
		double** chain1 = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			chain1[j] = new double[3];
			chain1[j][0] = coords[j + (i * chain_length)][0];
			chain1[j][1] = coords[j + (i * chain_length)][1];
			chain1[j][2] = coords[j + (i * chain_length)][2];
		}
		/**
		 *  Creating the entire chain2
		 */
		for (int j = i + 1; j < num_chains; j++)
		{
			double** chain2 = new double*[chain_length];
			/**
			 *  Scanning along both chains
			 */
			for (int k = 0; k < chain_length; k++)
			{
				chain2[k] = new double[3];
				chain2[k][0] = coords[k + (j * chain_length)][0];
				chain2[k][1] = coords[k + (j * chain_length)][1];
				chain2[k][2] = coords[k + (j * chain_length)][2];
			}

			/**
			 *  a is the scanning length which varies from start_chunk to end_chunk with a step size of
			 * step
			 */
			for (int a = start_chunk; a <= end_chunk; a += step)
			{
				double** temp_chain1 = new double*[a];
				for (int b = 0; b < a; b++)
				{
					temp_chain1[b] = new double[3];
				}
				for (int b = 0; b < chain_length - a; b++)
				{
					/**
					 *  Scanning along chain1
					 */
					for (int c = 0; c < a; c++)
					{
						temp_chain1[c][0] = chain1[b + c][0];
						temp_chain1[c][1] = chain1[b + c][1];
						temp_chain1[c][2] = chain1[b + c][2];
					}

					double** temp_chain2 = new double*[a];
					for (int e = 0; e < a; e++)
					{
						temp_chain2[e] = new double[3];
					}
					for (int e = 0; e < chain_length - a; e++)
					{
						/**
						 *  Scanning along chain2
						 */
						for (int f = 0; f < a; f++)
						{
							temp_chain2[f][0] = chain2[e + f][0];
							temp_chain2[f][1] = chain2[e + f][1];
							temp_chain2[f][2] = chain2[e + f][2];
						}

						/**
						 *  Calculating the LK between the two scanning instances
						 */
						double res = lk(temp_chain1, temp_chain2, a, a, false);
						outfile << "linking number between chains " << i << " starting at atom " << b << " and " << j
								<< " starting at atom " << e << " with chunk length " << a << ": " << res << "\n";
					}

					delete_array(temp_chain2, a);
				}

				delete_array(temp_chain1, a);
			}

			delete_array(chain2, chain_length);
		}

		delete_array(chain1, chain_length);
	}

	delete_array(coords, num_chains);
	outfile.close();

	return 0;
}
