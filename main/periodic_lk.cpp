#include "../include/funcs.h"

using namespace std;

int main(int argc, char* argv[])
{
	int num_chains = stoi(argv[3]);
	int chain_length = stoi(argv[2]);
	double box_dim;
	if (argc >= 5)
		box_dim = stod(argv[4]);
	else
		box_dim = 0;
	vector<double> box_dims = {box_dim, box_dim, box_dim};
	int num;
	double** coords = read_coords(argv[1], &num);
	double result[num_chains][num_chains];
	create_output_dir();
	ofstream outfile;
	outfile.open("./output/periodic_lk_out.txt");

	for (int i = 0; i < num_chains - 1; i++)
	{
		result[i][i] = 0;
		double** chain1 = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			chain1[j] = new double[3];
			chain1[j][0] = coords[j + (i * chain_length)][0];
			chain1[j][1] = coords[j + (i * chain_length)][1];
			chain1[j][2] = coords[j + (i * chain_length)][2];
		}
		// Images1 is the cells that chain1 intersects
		vector<vector<int>> images1 = compute_img(chain1, chain_length, box_dims);

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

			// Images2 is the cells that chain2 intersects
			vector<vector<int>> images2 = compute_img(chain2, chain_length, box_dims);
			vector<vector<int>> trans;
			for (int x = 0; x < images1.size(); x++)
			{
				for (int y = 0; y < images2.size(); y++)
				{
					vector<int> temp
						= {images1[x][0] - images2[y][0], images1[x][1] - images2[y][1], images1[x][2] - images2[y][2]};
					bool found = false;
					for (auto vec: trans)
					{
						if (vec == temp)
						{
							found = true;
							break;
						}
					}

					// trans contains all the images of chain2 that intersect the minmal unfolding of chain1
					// temp is one such image
					if (!found)
					{
						trans.push_back(temp);
						double res = lk(
							chain1,
							chain2,
							chain_length,
							chain_length,
							false,
							temp[0] * box_dims[0],
							temp[1] * box_dims[1],
							temp[2] * box_dims[2]);
						// at the end of the loop result[i][j] is the local periodic linking number between
						// chain i and chain j
						result[i][j] += res;
					}
				}
			}

			outfile << result[i][j] << "\n";
			delete_array(chain2, chain_length);
		}

		delete_array(chain1, chain_length);
	}

	delete_array(coords, num_chains);
	outfile.close();

	return 0;
}
