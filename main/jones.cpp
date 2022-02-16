#include "../include/funcs.h"
using namespace std;

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
		int loops = count_loops(neigh_array, n);
		if (loops == 1)
		{
			map<int, double> res;
			res[0] = 1;
			return res;
		}
		else
		{
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

		if (compute_one(coords[p1], coords[p2], coords[p3], coords[p4]) < 0)
		{
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

		return res;
	}
}

map<int, double> jones(string fname, double** coordinates, bool has_coords, int length, int num_proj)
{
	int n = length;
	bool is_closed = false;
	vector<vector<double>> coords;
	for (int i = 0; i < n; i++)
	{
		coords.push_back({coordinates[i][0], coordinates[i][1], coordinates[i][2]});
	}
	n = coords.size();

	vector<vector<int>> neigh_array = generate_neigh_array(n, is_closed);
	map<int, double> avejones;
	vector<vector<double>> original_coords(coords);
	int num_fails = 0;
	for (int z = 0; z < num_proj; z++)
	{
		coords = original_coords;
		n = coords.size();
		neigh_array = generate_neigh_array(n, is_closed);
		bool no_mult_crossings = false;
		vector<vector<double>> proj;
		vector<vector<double>> new_proj = get_proj(coords, n);
		vector<vector<int>> crossings = count_crossings(neigh_array, new_proj, n, false);
		Struct s = next_mult_crossings(new_proj, coords, neigh_array, n, is_closed);
		if (!s.success)
		{
			z--;
			continue;
		}
		if (s.crossings.size() > crossings.size() || has_mult_crossings(s.crossings))
		{
			z--;
			continue;
		}

		proj = s.proj;
		coords = s.coords;
		crossings = s.crossings;
		n = coords.size();
		neigh_array = generate_neigh_array(n, is_closed);
		bool many_crossings = false;
		int init_crossings = crossings.size();

		vector<vector<int>> temp_crossings = reduce_crossings(crossings, coords);
		crossings = temp_crossings;
		init_crossings = temp_crossings.size();
		if (init_crossings > 15)
		{
			z--;
			num_fails++;
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
				continue;
			}
		}

		int writhe = 0;
		for (int i = 0; i < init_crossings; i++)
		{
			if (are_collinear(
					coords[crossings[i][0]], coords[crossings[i][1]], coords[crossings[i][2]], coords[crossings[i][3]]))
				continue;
			if (compute_one(
					coords[crossings[i][0]], coords[crossings[i][1]], coords[crossings[i][2]], coords[crossings[i][3]])
				> 0)
				writhe++;
			else
				writhe--;
		}

		map<int, double> result = anneal(neigh_array, coords, proj, crossings, init_crossings, n);
		map<int, double> mult;
		mult[(-1 * writhe) * 3] = pow(-1, (-1 * writhe) * 3);
		map<int, double> mult_result = simple_mult(result, mult);
		for (map<int, double>::const_iterator it = mult_result.begin(); it != mult_result.end(); ++it)
		{
			if (it->second == 0)
				continue;
			avejones[it->first] += it->second;
		}
	}

	map<int, double> final;
	for (map<int, double>::const_iterator it = avejones.begin(); it != avejones.end(); ++it)
	{
		final[it->first] = it->second / num_proj;
	}

	return final;
}

int main(int argc, char* argv[])
{
	int num_chains = stoi(argv[3]);
	int chain_length = stoi(argv[2]);
	// int file_type = stoi(argv[4]);
	int num_repeats = 5;
	int num;
	double box_dim;
	if (argc >= 5)
		box_dim = stod(argv[4]);
	else
		box_dim = 0;
	vector<double> box_dims = {box_dim, box_dim, box_dim};
	ofstream outfile;
	create_output_dir();
	outfile.open("./output/jones_out.txt");
	map<int, double> result;
	double** coords = read_coords(argv[1], &num);
	for (int i = 0; i < num_chains; i++)
	{
		double** temp_coords = new double*[chain_length];
		for (int j = 0; j < chain_length; j++)
		{
			temp_coords[j] = new double[3];
			temp_coords[j][0] = coords[(i * chain_length) + j][0];
			temp_coords[j][1] = coords[(i * chain_length) + j][1];
			temp_coords[j][2] = coords[(i * chain_length) + j][2];
		}

		map<int, double> jones_poly = jones("", temp_coords, true, chain_length, 10);
		for (map<int, double>::const_iterator it = jones_poly.begin(); it != jones_poly.end(); ++it)
		{
			outfile << it->second << "A^" << it->first << " + ";;
			result[it->first] += it->second;
		}
		outfile << "\n";
		delete_array(temp_coords, chain_length);
	}

	for (map<int, double>::const_iterator it = result.begin(); it != result.end(); ++it)
	{
		if (it->second == 0 || it->second == -0)
			continue;
		outfile << it->second / num_chains << "A^" << it->first << " + ";
	}
	outfile << "\n";

	outfile.close();
	return 0;
}
