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
using namespace std;

/* Given the chain coordinates teh periodic box size and whether is it open or closed, returns the periodic writhe
*/
double periodic_wr(double** coords, int chain_length, vector<double> box_dims, bool is_closed)
{
	double result = 0;
	vector<vector<int>> images = compute_img(coords, chain_length, box_dims);
        vector<vector<int>> allimages = compute_periodic_img(images, coords, chain_length, box_dims);

	for (int i = 0; i < allimages.size(); i++)
	{
  
	    if (allimages[i][0] == 0 && allimages[i][1] == 0 && allimages[i][2] == 0)
			continue;
            else
                {
                  vector<int> translation_vector = {allimages[i][0], allimages[i][1], allimages[i][2]};
                  //vector<vector<double>> 
                  double** coords2 = translation(coords, chain_length,  translation_vector, box_dims);
                  result += lk(coords,coords2, chain_length,chain_length,is_closed,0,0,0); 
                }
               

	}

	return result;
}

/*input: the filename, the chain length, the number of chains, 0/1 for ring/linear, box dimensions
returns the periodic writhe of all chains in the system
*/
int main(int argc, char* argv[])
{
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        //cout << "ENTER";
	int num_chains = stoi(argv[3]);
	int chain_length = stoi(argv[2]);
        int chunk = num_chains / size;
	int num;
        int ringlinear=stoi(argv[4]);       
        bool is_closed;
	double box_dim;
	if (argc >= 6)
		box_dim = stod(argv[5]);
	else
		box_dim = 0;
        //cout<< box_dim;
        if (ringlinear==0)
        {
           is_closed=true;
        }
        else
        {
           is_closed=false;
        }

	vector<double> box_dims = {box_dim, box_dim, box_dim};
        create_output_dir();
        if (!fs::exists("./output/periodic_wr_mpi"))
        {
                cout << "Creating periodic_wr_mpi directory..." << endl;
                fs::create_directory("./output/periodic_wr_mpi");
        }
        ofstream outfile;
        string file_name = to_string(chain_length) + "_periodic_wr_mpi_out_" + to_string(rank) + ".txt";
        outfile.open("./output/periodic_wr_mpi/" + file_name);
	//create_output_dir();
	//ofstream outfile;
	//outfile.open("./output/periodic_wr_out.txt");
	double** coords = read_coords(argv[1], &num, chain_length, box_dim);
	for (int i = rank*chunk; i < (rank+1)*chunk; i++)
	{
		double** temp_coords = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			temp_coords[j] = new double[3];
			temp_coords[j][0] = coords[(i * chain_length) + j][0];
			temp_coords[j][1] = coords[(i * chain_length) + j][1];
			temp_coords[j][2] = coords[(i * chain_length) + j][2];
			// cout << temp_coords[j][0] << ", " << temp_coords[j][1] << ", " << temp_coords[j][2] <<
			// "\n";
		}

		double pwr = periodic_wr(temp_coords, chain_length, box_dims, is_closed);
		//pwr += wr(temp_coords, chain_length, is_closed);
		outfile << pwr << "\n";
		delete_array(temp_coords, chain_length);
	}

	outfile.close();
	delete_array(coords, num_chains);
        MPI_Finalize();
	return 0;
}
