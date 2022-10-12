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

/*Given a projection of a chain and its crossings, anneal crossings using the Skein Relations
*/
map<int, double> anneal(
	vector<vector<int>> neigh_array,
	vector<vector<double>> coords,
	vector<vector<double>> proj,
	vector<vector<int>> crossings,
	int init_crossings,
	int n)
{
	if (crossings.size() == 0)
	{
                //cout<<"No crossings left \n";
		int loops = count_loops(neigh_array, n);
		if (loops == 1)
		{
                        //cout<<"ONLY ONE LOOP \n";
			map<int, double> res;
			res[0] = 1;
			return res;
		}
		else
		{
                        //cout<<"LOOPS = "<<loops<<"\n";
			map<int, double> res = mult_poly(loops);
			if (res.find(-5) != res.end())
				cout << "0 crossings left\n";
			return res;
		}
	}
	else
	{
		int p1 = crossings[0][0];
		int p2 = crossings[0][1];
		int p3 = crossings[0][2];
		int p4 = crossings[0][3];
		vector<vector<int>> neigh_copy1(neigh_array);
		vector<vector<int>> neigh_copy2(neigh_array);
                //cout << "We are annealing crossings "<< p1<<","<<p2<<","<<p3<<","<<p4<<"\n";
		if (compute_one(coords[p1], coords[p2], coords[p3], coords[p4]) < 0)
		{
                        //cout<<"They have negative Writhe \n";
			neigh_copy1[p1][1] = neigh_array[p3][1];
			neigh_copy1[p3][1] = neigh_array[p1][1];
			neigh_copy1[neigh_array[p1][1]][0] = p3;
			neigh_copy1[neigh_array[p3][1]][0] = p1;
			neigh_copy2[p1][1] = p3;
			neigh_copy2[p3][1] = p1;
			neigh_copy2[neigh_array[p1][1]][0] = neigh_array[p3][1];
			neigh_copy2[neigh_array[p3][1]][0] = neigh_array[p1][1];
		}
		else
		{
                        //cout<<"They have positive writhe \n";
			neigh_copy2[p1][1] = neigh_array[p3][1];
			neigh_copy2[p3][1] = neigh_array[p1][1];
			neigh_copy2[neigh_array[p1][1]][0] = p3;
			neigh_copy2[neigh_array[p3][1]][0] = p1;
			neigh_copy1[p1][1] = p3;
			neigh_copy1[p3][1] = p1;
			neigh_copy1[neigh_array[p1][1]][0] = neigh_array[p3][1];
			neigh_copy1[neigh_array[p3][1]][0] = neigh_array[p1][1];
		}

		crossings.erase(crossings.begin());
		map<int, double> ret1 = anneal(neigh_copy1, coords, proj, crossings, init_crossings, n);
		map<int, double> ret2 = anneal(neigh_copy2, coords, proj, crossings, init_crossings, n);
		map<int, double> temp1;
		map<int, double> temp2;
		temp1[-1] = 1;
		temp2[1] = 1;
		map<int, double> res1 = simple_mult(ret1, temp1);
		map<int, double> res2 = simple_mult(ret2, temp2);
		map<int, double> res;
		for (auto const& x: res1)
		{
			res[x.first] += x.second;
		}
		for (auto const& x: res2)
		{
			res[x.first] += x.second;
		}
                //cout<<"The result so far is ",res[0],",",res[1],"\n";
		return res;
	}
}

/* Given the chain coordinates, whether it is closed or open and the number of projections,
compute the Jones polynomial of the chain
*/
map<int, double> jones(string fname, double** coordinates, int length, bool is_closed, int num_proj)
{
	int n = length;
	vector<vector<double>> coords;
	for (int i = 0; i < n; i++)
	{
		coords.push_back({coordinates[i][0], coordinates[i][1], coordinates[i][2]});
	}


	vector<vector<int>> neigh_array = generate_neigh_array(n, is_closed);
	map<int, double> avejones;
	int num_fails = 0;
        //double twr=wr(coordinates,n,is_closed);
        if (is_closed)
        {
            int m = 4;
        }
        //int exclude = 0;
        int init_n=n;
        int tot=0;
        vector<vector<double>> init_coords=coords;
        vector<vector<int>> init_neigh_array=neigh_array;
	for (int z = 0; z < num_proj; z++)
	{
                bool multc = true;
                vector<vector<double>> proj;
                vector<vector<double>> new_proj = get_proj(init_coords,init_n);
                vector<vector<int>> crossings = count_crossings(init_neigh_array, new_proj, init_n, is_closed);
                if (crossings.size()==0)
                       continue;
               
                  if (multc==true)
                  {
                  
                    int initcross=crossings.size();
                    Struct s = mult_crossings(new_proj, init_coords, init_neigh_array, crossings, init_n, is_closed);
		    if (!s.success)
		    {
			z--;
                        //exclude++;
			continue;
		    }

		    new_proj = s.proj;
		    coords = s.coords;
                
	            n = coords.size();
                    
		    neigh_array = generate_neigh_array(n, is_closed);
                    crossings = count_crossings(neigh_array, new_proj, n, is_closed);
                    
                    int aftercross=crossings.size();
                    if (initcross!=aftercross)
                    {
                       z--;
                       //exclude++;
                       continue;
                    }
                    multc=mult_cross(n,crossings,is_closed);
                
                }
		bool many_crossings = false;
		int init_crossings = crossings.size();
                
		vector<vector<int>> temp_crossings = reduce_crossings(crossings, coords);
		crossings = temp_crossings;
                
		init_crossings = temp_crossings.size();
		if (init_crossings > 15)
		{
			z--;
			num_fails++;
                        //exclude++;
			if (num_fails >= 25)
			{
				map<int, double> err_map;
				err_map[0] = 0;
                                
				return err_map;
			};
			continue;
			many_crossings = true;

			if (init_crossings > 15)
			{
				z--;
                                //exclude++;
				continue;
			}
		}

		int writhe = 0;
		for (int i = 0; i < init_crossings; i++)
		{
			
			if (compute_one(
					coords[crossings[i][0]], coords[crossings[i][1]], coords[crossings[i][2]], coords[crossings[i][3]])
				> 0)
				writhe++;
			else
				writhe--;
		}
                

		map<int, double> result = anneal(neigh_array, coords, new_proj, crossings, init_crossings, n);
		map<int, double> mult;
		mult[(-1 * writhe) * 3] = pow(-1, (-1 * writhe));
		map<int, double> mult_result = simple_mult(result, mult);
                tot++;
		for (map<int, double>::const_iterator it = mult_result.begin(); it != mult_result.end(); ++it)
		{
			if (it->second == 0)
                                
				continue;
			avejones[it->first] += it->second;
		}
                proj.clear();
                new_proj.clear();
                coords.clear();
                neigh_array.clear();
	}

	map<int, double> final;
	for (map<int, double>::const_iterator it = avejones.begin(); it != avejones.end(); ++it)
	{
                
		final[it->first] = it->second / tot;
	}

	return final;
}

/* -*- -*- ----------------------------------------------------------
   Takes as input the filename, the number of chains, the length of the chains, integer 0/1 statement (depending on whether the chains are closed, i.e. rings, or open, i.e. linear, 0: ring, 1: linear), the numbber of projections,the start chunk, the end chunk, the chunk step,  and the box dimension (optional)
   If the box dimension is not specified, or if it is equal to 0, then the system is not periodic and the coordinates are unwrapped. 
   If a non-zero box-dimension is specified, the coordinates are unwrapped, according to the PBC.
   Returns the Jones polynomial of each subsegement of chunk size along each chain
------------------------------------------------------------------------- */

int main(int argc, char* argv[])
{
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int num_chains = stoi(argv[3]);
	int chain_length = stoi(argv[2]);
	
       
        bool is_closed;
        int chunk = num_chains / size;
	int num;
        int ringlinear =stoi(argv[4]);
        int nproj = stoi(argv[5]);
        int start_chunk=stoi(argv[6]);
        int end_chunk=stoi(argv[7]);
        int step=stoi(argv[8]);
	double box_dim;
	if (argc >= 10)
        {
		box_dim = stod(argv[9]);
        }
	else
        {
		box_dim = 0;
        }
	vector<double> box_dims = {box_dim, box_dim, box_dim};
	//ofstream outfile;
        if (ringlinear==0)
        {
           is_closed = true;
        }
        else
        {
           is_closed = false;
        }
        create_output_dir();
        string file_name = to_string(chain_length) + "_jones_scan_mpi_out_" + to_string(rank) + ".txt";
        if (!fs::exists("./output/jones_scan_mpi"))
        {
                cout << "Creating jones_scan_mpi directory..." << endl;
                fs::create_directory("./output/jones_scan_mpi");
        }
        ofstream outfile;
        outfile.open("./output/jones_scan_mpi/" + file_name);
	//create_output_dir();
	//outfile.open("./output/jones_scan_out.txt");
	//map<int, double> result;
	double** coords; //= read_coords(argv[1], &num);
        if (box_dim == 0)
        {
                coords = read_coords(argv[1], &num);
        }
        else
        {
                coords = read_coords(argv[1], &num, chain_length, box_dim);
        }
	for (int i = rank*chunk; i < (rank+1)*chunk; i++)
	{      
           
                //cout<< "COMPUTING JONES OF CHAIN "<<i<<"------------------------------------------------------------- \n";
		double** temp_coords = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			temp_coords[j] = new double[3];
			temp_coords[j][0] = coords[(i * chain_length) + j][0];
			temp_coords[j][1] = coords[(i * chain_length) + j][1];
			temp_coords[j][2] = coords[(i * chain_length) + j][2];
		}
                for (int a=start_chunk; a<=end_chunk; a+=step)
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
                                        temp_chain1[c][0] = temp_coords[b + c][0];
                                        temp_chain1[c][1] = temp_coords[b + c][1];
                                        temp_chain1[c][2] = temp_coords[b + c][2];
                                }
                                map<int, double> result;
		                map<int, double> jones_poly = jones("", temp_chain1, a,is_closed,  nproj);
		                for (map<int, double>::const_iterator it = jones_poly.begin(); it != jones_poly.end(); ++it)
	                	{
			            outfile << it->second << "A^" << it->first << " + ";;
			            //result[it->first] += it->second;
	                	}
		                outfile << "\n";
		                //delete_array(temp_chain1, a);
                    }
                    delete_array(temp_chain1,a);
                 }
                 delete_array(temp_coords, chain_length);

	}

	//for (map<int, double>::const_iterator it = result.begin(); it != result.end(); ++it)
	//{
	//	if (it->second == 0 || it->second == -0)
	//		continue;
        //	outfile << it->second / num_chains << "A^" << it->first << " + ";
	//}
        delete_array(coords, num_chains);
	outfile << "\n";

	outfile.close();
        MPI_Finalize();
	return 0;
}
