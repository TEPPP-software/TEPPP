#ifndef FUNCS_H
#define FUNCS_H

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>
namespace fs = std::filesystem;
using namespace std;




//random number generator


#ifdef __cplusplus
extern "C"
{
#endif
#define CHECK(x)                \
	do                          \
	{                           \
		if (0 != (error = (x))) \
			goto fail;          \
	} while (0)
#define FAIL(...)                     \
	do                                \
	{                                 \
		fprintf(stderr, __VA_ARGS__); \
		error = 1;                    \
		goto fail;                    \
	} while (0)

	typedef struct
	{
		FILE* file;
		uint64_t* array;
		int offset;
	} ruv_buffer_t;

	int ruv_init(ruv_buffer_t* buffer);
	int ruv_uniform(ruv_buffer_t* buffer, double* x);
	int ruv_generate(ruv_buffer_t* buffer, int ndimensions, double vector[]);
	void ruv_free(ruv_buffer_t* buffer);

	static const double TWO_PI = 6.283185307179586476925287;

#define BUFFER_LENGTH 1024

	static const char* RandomFileName = "/dev/urandom";

	int ruv_init(ruv_buffer_t* buffer)
	{
		int error;
		memset(buffer, 0, sizeof(ruv_buffer_t));

		buffer->file = fopen(RandomFileName, "rb");
		if (buffer->file == NULL)
			FAIL("ruv_init: cannot open random device %s\n", RandomFileName);

		buffer->array = static_cast<uint64_t*>(calloc(BUFFER_LENGTH, sizeof(buffer->array[0])));
		if (buffer->array == NULL)
			FAIL("ruv_init: cannot allocate memory\n");

		/* Cause first call to ruv_uniform() to trigger a read. */
		buffer->offset = BUFFER_LENGTH;

		error = 0;
	fail:
		if (error != 0)
			ruv_free(buffer);

		return error;
	}


	int ruv_uniform(ruv_buffer_t* buffer, double* x)
	{
		int error;

		/* Keep reading until the resulting random number is not zero. */
		do
		{
			if (buffer->offset == BUFFER_LENGTH)
			{
				/* We need to load a new batch of random 64-bit unsigned integers. */
				size_t nread = fread(buffer->array, sizeof(buffer->array[0]), BUFFER_LENGTH, buffer->file);
				if (nread != BUFFER_LENGTH)
					FAIL("ruv_uniform: Error reading from device %s\n", RandomFileName);
				buffer->offset = 0;
			}

			*x = (double) buffer->array[buffer->offset++] / (double) 0xffffffffffffffffLU;
		} while (*x == 0.0);

		/* Now *x contains a uniformly-distributed random value in the half-open range (0, 1]. */
		error = 0;
	fail:
		return error;
	}


	int ruv_generate(ruv_buffer_t* buffer, int ndimensions, double vector[])
	{
		int error, i;
		double A, B; /* uniform random variables */
		double radius, angle;
		double sum;

		do
		{
			sum = 0.0;
			for (i = 0; i < ndimensions; i += 2)
			{
				CHECK(ruv_uniform(buffer, &A));
				CHECK(ruv_uniform(buffer, &B));
				radius = sqrt(-2 * log(A));
				angle = TWO_PI * B;
				vector[i] = radius * cos(angle);
				sum += vector[i] * vector[i];
				if (i + 1 < ndimensions)
				{
					vector[i + 1] = radius * sin(angle);
					sum += vector[i + 1] * vector[i + 1];
				}
			}
		} while (sum == 0.0);

		/* Convert to a unit vector by dividing through by the length. */
		sum = sqrt(sum);
		for (i = 0; i < ndimensions; ++i)
			vector[i] /= sum;

		error = 0;
	fail:
		return error;
	}


	void ruv_free(ruv_buffer_t* buffer)
	{
		if (buffer->file != NULL)
		{
			fclose(buffer->file);
			buffer->file = NULL;
		}

		if (buffer->array != NULL)
		{
			free(buffer->array);
			buffer->array = NULL;
		}
	}

#ifdef __cplusplus
}
#endif

struct proj_coords
{
	vector<vector<double>> coords;
	vector<vector<double>> proj;
	vector<vector<int>> crossings;
	bool success;
};

typedef struct proj_coords Struct;

void delete_array(double** arr, int size);
void delete_array(int** arr, int size);
void normalize3(const double* v, double* ans);
void cross3(const double* v1, const double* v2, double* ans);
void sub3(const double* v1, const double* v2, double* ans);

bool intersect(double* p1, double* q1, double* p2, double* q2);
bool are_collinear(vector<double> p1, vector<double> p2, vector<double> p3, vector<double> p4);
bool intersect1(vector<double> p0, vector<double> p1, vector<double> p2, vector<double> p3, double* rx, double* ry);
map<int,vector<int>> has_mult_crossings(int n, vector<vector<int>> crossings, bool is_closed);

int count_loops(vector<vector<int>> neigh_array, int n);

double distance(vector<double> p1, vector<double> p2);
double compute_one(vector<double> p1, vector<double> p2, vector<double> p3, vector<double> p4);
double dot3(const double* v1, const double* v2);
double wr(double** chain, int length, bool is_closed);
double
lk(double** chain1,
   double** chain2,
   int length1,
   int length2,
   bool is_closed,
   double offx = 0,
   double offy = 0,
   double offz = 0);

vector<vector<double>> get_random_proj();
vector<double> unwrap(vector<double> p1, vector<double> p2, double box_dim);
vector<double> get_intersection(vector<double> p0, vector<double> p1, vector<double> p2, vector<double> p3);
map<int, double> mult_poly(int power);
map<int, double> simple_mult(map<int, double> a, map<int, double> b);

vector<vector<int>> generate_neigh_array(int n, bool closed);
vector<vector<int>> compute_img(double** coords, int n, vector<double> box_dims);
vector<vector<int>> compute_periodic_img(vector<vector<int>> cells, double** coords, int n, vector<double> box_dims);
//vector<vector<double>> 
double** translation(double** coords,int n,  vector<int> translation_vector, vector<double> box_dims );
vector<vector<int>>
count_crossings(vector<vector<int>> neigh_array, vector<vector<double>> proj, int n, bool count_all);
vector<vector<int>> reduce_crossings(vector<vector<int>> initial_crossings, vector<vector<double>> coords);

double** read_coords(string fname, int* n);
double** read_coords(string fname, int* n, int chain_length, double box_dim);
vector<vector<double>> invert_mat(vector<vector<double>> m);
vector<vector<double>> get_proj(vector<vector<double>> coords, int n);

Struct next_mult_crossings(
	vector<vector<double>> proj, vector<vector<double>> coords, vector<vector<int>> neigh_array, int n, bool is_closed);

Struct mult_crossings(
        vector<vector<double>> proj, vector<vector<double>> coords, vector<vector<int>> neigh_array, vector<vector<int>> before_cross, int n, bool is_closed);


vector<int> find_knot(map<int, double> poly, vector<string> params);

/**
 * Delete 2D array of doubles
 *
 * @param arr Array of doubles to delete
 * @param size Number of elements in array
 *
 */
void delete_array(double** arr, int size)
{
	for (int i = 0; i < size; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;
}

/**
 * Delete 2D array of ints
 *
 * @param arr Array of ints to delete
 * @param size Number of elements in array
 *
 */
void delete_array(int** arr, int size)
{
	for (int i = 0; i < size; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;
}

/**
 * Normalize a 3D vector
 *
 * @param[in] v The vector to normalize
 * @param[out] ans The normalized vector
 */
void normalize3(const double* v, double* ans)
{
	double scale = 1.0 / sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	ans[0] = v[0] * scale;
	ans[1] = v[1] * scale;
	ans[2] = v[2] * scale;
}

/**
 * Get the cross product between two 3D vectors
 *
 * @param[in] v1 The first vector in the calculation
 * @param[in] v2 The second vector in the calculation
 * @param[out] ans The cross product of v1 and v2
 */
void cross3(const double* v1, const double* v2, double* ans)
{
	ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
	ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
	ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/**
 * Get the difference between two 3D vectors
 *
 * @param[in] v1 The first vector in the calculation
 * @param[in] v2 The second vector in the calculation
 * @param[out] ans The difference of v1 and v2
 */
void sub3(const double* v1, const double* v2, double* ans)
{
	ans[0] = v1[0] - v2[0];
	ans[1] = v1[1] - v2[1];
	ans[2] = v1[2] - v2[2];
}

/**
 * Check if 4 3D points are collinear
 *
 * @param p1, p2, p3, p4 Points to check for collinearity
 * @return true if points are collinear, false if points are not collinear
 */
bool are_collinear(vector<double> p1, vector<double> p2, vector<double> p3, vector<double> p4)
{
	vector<double> v1 = {p2[0] - p1[0], p2[1] - p1[1]};
	vector<double> v2 = {p4[0] - p3[0], p4[1] - p3[1]};
	double cross = v1[0] * v2[1] - v2[0] * v1[1];

	if (abs(cross) < 0.001)
		return true;
	else
		return false;
}

/**
 * Check if two edges in 3-space intersect
 *
 * @param[in] p0, p1 Points that form the first edge
 * @param[in] p2, p3 Points that form the second edge
 * @param[out] rx, ry Coordinates of intersection point if the edges intersect
 * @return true if edges intersect, false if edges do not intersect
 */
bool intersect1(vector<double> p0, vector<double> p1, vector<double> p2, vector<double> p3, double* rx, double* ry)
{
	if (are_collinear(p0, p1, p2, p3))
		return false;
	double s1_x, s1_y, s2_x, s2_y;
	s1_x = p1[0] - p0[0];
	s1_y = p1[1] - p0[1];
	s2_x = p3[0] - p2[0];
	s2_y = p3[1] - p2[1];

	double det = -s2_x * s1_y + s1_x * s2_y;
	if (abs(det) == 0.0)
		return false;
	double s, t;
	s = (-s1_y * (p0[0] - p2[0]) + s1_x * (p0[1] - p2[1])) / det;
	t = (s2_x * (p0[1] - p2[1]) - s2_y * (p0[0] - p2[0])) / det;

	if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
	{
		if (rx != NULL)
		{
			*rx = p0[0] + (t * s1_x);
		}

		if (ry != NULL)
		{
			*ry = p0[1] + (t * s1_y);
		}

		return true;
	}

	return false;
}

/**
 * Creates a list of edges that each edge intersects
 *
 * @param crossings Vector of vectors each containing the indices of four points that form two crossing edges
 * @return mult: the list of edges that each edge intersects
 */
map<int, vector<int>>  has_mult_crossings(int n,vector<vector<int>> crossings, bool is_closed)
{
	map<int, bool> found;
        map<int, bool> found2;
        map<int, vector<int>> mult;
        int last = n;
        if (is_closed)
        {
           last=n+1;
        }
        else
        {
           last=n;
        }
        for (int i=0;i<last;i++)
        {
                mult[i]={i};
                found[i]=false;
                found2[i]=false;       
        }
        

        vector<vector<int>> cross3;
	for (int i = 0; i < last; i++) 
	{
             
                vector<int> cross;
                vector<int> cross2; 
                int a = i;
                vector<int> crosslist;
                for (int j=0; j<crossings.size();j++)
                {
                      int c=crossings[j][0];
                      int d=crossings[j][2];
                      if (a==c)
                         
                         cross.push_back(d);
                                        
                      if (a==d)
                         
                         cross.push_back(c);
                                 
                }
                mult[i]=cross;
                
                
	}
  
        
        for(auto itr = mult.begin(); itr != mult.end(); itr++) 
        {
             int i = itr->first;
             //cout<<"length "<< mult[i].size()<<"\n";
             
        }
        cout << endl;

	return mult;
}


/**
 * Checks if there is an edge that intersects with more than one other edge
 *
 * @param mult contains the list of edges each edge intersect
 * @return true if there is an edge that intesects more than one other edge, false if not
 */


bool mult_cross(int n,vector<vector<int>> crossings, bool is_closed)
{
  map<int, vector<int>> mult=has_mult_crossings(n,crossings,is_closed);
  bool res=false;
  for (int i=0;i<mult.size();i++)
  {
     if (mult[i].size()>1)
        res=true;
  }

  return res;
}
/**
 * Counts the number of loops formed by order of atoms in neigh_array.
 * Used to calculate the bracket polynomial of each combination of crossing annealments.
 *
 * @param neigh_array Vector of vectors that contain an the indices of an atom's neighbors
 * @param n The number of atoms in the system
 * @return The number of loops formed by the system
 */
int count_loops(vector<vector<int>> neigh_array, int n)
{
	bool found[n] = {false};
	int count = 0;
	int loops = 0;
	for (int i = 0; i < n; i++)
	{
		if (found[i])
			continue;
		found[i] = true;
		int curr = i;
		while (1)
		{
			found[curr] = true;
			if (neigh_array[curr][0] >= 0 && !found[neigh_array[curr][0]])
			{
				int temp = neigh_array[curr][0];
				curr = temp;
			}
			else if (neigh_array[curr][1] >= 0 && !found[neigh_array[curr][1]])
			{
				int temp = neigh_array[curr][1];
				curr = temp;
			}
			else
			{
				loops++;
				break;
			}
		}
	}

	return loops;
}

/**
 * Gives the distance between two points
 *
 * @param p1, p2 The points to measure the distance between
 * @return The distance between p1 and p2
 */
double distance(vector<double> p1, vector<double> p2) { return sqrt(pow(p2[0] - p1[0], 2) + pow(p2[1] - p1[1], 2)); }

/**
 * Compute the linking number between two edges formed by two points each
 *
 * @param p1, p2 The two points that form the first edge in the calculation
 * @param p3, p4 The two points that form the second edge in the calculation
 * @return The linking number of the two edges formed by p1, p2, p3, and p4
 */
double compute_one(vector<double> p1, vector<double> p2, vector<double> p3, vector<double> p4)
{
	double a1[3];
	double a2[3];
	double b1[3];
	double b2[3];
	double c1[3];
	double sum;
	double ra[3];
	double rb[3];
	double r00[3];
	double r01[3];
	double r10[3];
	double r11[3];
	double v1[3];
	double v2[3];
	double v3[3];
	double v4[3];
	double u1[3];
	double u2[3];
	double u3[3];
	double u4[3];
	double d1;
	double d2;
	double d3;
	double d4;
	double as1;
	double as2;
	double as3;
	double as4;
	double aux1[3];
	double aux;
	double alk;
	double norm;
	double sign;
	double pi = 2 * asin(1.0);

	a1[0] = p1[0];
	a1[1] = p1[1];
	a1[2] = p1[2];
	a2[0] = p2[0];
	a2[1] = p2[1];
	a2[2] = p2[2];
	b1[0] = p3[0];
	b1[1] = p3[1];
	b1[2] = p3[2];
	b2[0] = p4[0];
	b2[1] = p4[1];
	b2[2] = p4[2];
	sub3(a2, a1, ra);
	sub3(b2, b1, rb);
	sub3(a1, b1, r00);
	sub3(a1, b2, r01);
	sub3(a2, b1, r10);
	sub3(a2, b2, r11);
	cross3(r00, r01, v1);
	normalize3(v1, u1);
	cross3(r01, r11, v2);
	normalize3(v2, u2);
	cross3(r11, r10, v3);
	normalize3(v3, u3);
	cross3(r10, r00, v4);
	normalize3(v4, u4);
	d1 = dot3(u1, u2);
	d2 = dot3(u2, u3);
	d3 = dot3(u3, u4);
	d4 = dot3(u4, u1);
	as1 = asin(d1);
	as2 = asin(d2);
	as3 = asin(d3);
	as4 = asin(d4);
	cross3(ra, rb, aux1);
	aux = dot3(aux1, r00);
	if (aux < 0)
	{
		alk = -1 * ((as1 + as2 + as3 + as4) / (4 * pi));
	}
	else
	{
		alk = (as1 + as2 + as3 + as4) / (4 * pi);
	}
	if (alk != alk)
		return 0;
	else
		return alk;
}

/**
 * Compute the dot product between two 3d vectors
 */
double dot3(const double* v1, const double* v2) { return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]; }

/**
 * Compute the writhe of a chain
 *
 * @param chain A 2d array containing the coordinates for the chain to compute the writhe of
 * @param length The number of atoms in the chain
 * @param is_closed True if the first atom of the chain is connected to the last, false otherwise
 * @return The writhe of the chain
 */
double wr(double** chain, int length, bool is_closed)
{
	double result = 0;
	double pi = 2 * asin(1.0);

	for (int i = 0; i < length - 2; i++)
	{
		vector<double> p1, p2;
		if (i + 1 == length)
		{
			if (!is_closed)
				continue;
			p1.push_back(chain[i][0]);
			p1.push_back(chain[i][1]);
			p1.push_back(chain[i][2]);
			p2.push_back(chain[0][0]);
			p2.push_back(chain[0][1]);
			p2.push_back(chain[0][2]);
			// continue;
		}
		else
		{
			p1.push_back(chain[i][0]);
			p1.push_back(chain[i][1]);
			p1.push_back(chain[i][2]);
			p2.push_back(chain[i + 1][0]);
			p2.push_back(chain[i + 1][1]);
			p2.push_back(chain[i + 1][2]);
		}
		for (int j = i + 2; j < length; j++)
		{
			vector<double> p3, p4;
			if (j + 1 == length)
			{
				if (!is_closed)
					continue;
				p3.push_back(chain[j][0]);
				p3.push_back(chain[j][1]);
				p3.push_back(chain[j][2]);
				p4.push_back(chain[0][0]);
				p4.push_back(chain[0][1]);
				p4.push_back(chain[0][2]);
				// continue;
			}

			else
			{
				p3.push_back(chain[j][0]);
				p3.push_back(chain[j][1]);
				p3.push_back(chain[j][2]);
				p4.push_back(chain[j + 1][0]);
				p4.push_back(chain[j + 1][1]);
				p4.push_back(chain[j + 1][2]);
			}
                        double lwr=compute_one(p1,p2,p3,p4);
                        //cout<<"EDGE "<<i<<" and EDGE "<< j<<"WRITHE = "<<lwr<<"\n";
			result += compute_one(p1, p2, p3, p4);
		}
	}

	return result * 2;
}

/**
 * Compute the linking number between two chains
 *
 * @param chain1 The coordinates of the first chain in the calculation
 * @param chain2 The coordinates of the second chain in the calculation
 * @param length1 The number of atoms in chain1
 * @param length2 The number of atoms in chain2
 * @param is_closed True if the last atom of the chains are connected to the first, false otherwise
 * @param offx, offy, offz The amount of offset applied to the coordinates of the second chain (used in periodic linking
 * number calculation)
 * @return The linking number between chain1 and chain2
 */
double
lk(double** chain1,
   double** chain2,
   int length1,
   int length2,
   bool is_closed,
   double offx /*= 0*/,
   double offy /*= 0*/,
   double offz /*= 0*/)
{
	vector<double> p1, p2, p3, p4;
	double result = 0;
	double pi = 2 * asin(1.0);

	for (int i = 0; i < length1; i++)
	{
		vector<double> p1, p2;
		if (i + 1 == length1)
		{
			if (!is_closed)
				continue;
			p1.push_back(chain1[i][0]);
			p1.push_back(chain1[i][1]);
			p1.push_back(chain1[i][2]);
			p2.push_back(chain1[0][0]);
			p2.push_back(chain1[0][1]);
			p2.push_back(chain1[0][2]);
			// continue;
		}
		else
		{
			p1.push_back(chain1[i][0]);
			p1.push_back(chain1[i][1]);
			p1.push_back(chain1[i][2]);
			p2.push_back(chain1[i + 1][0]);
			p2.push_back(chain1[i + 1][1]);
			p2.push_back(chain1[i + 1][2]);
		}
		for (int j = 0; j < length2; j++)
		{
			vector<double> p3, p4;
			if (j + 1 == length2)
			{
				if (!is_closed)
					continue;
				p3.push_back(chain2[j][0] + offx);
				p3.push_back(chain2[j][1] + offy);
				p3.push_back(chain2[j][2] + offz);
				p4.push_back(chain2[0][0] + offx);
				p4.push_back(chain2[0][1] + offy);
				p4.push_back(chain2[0][2] + offz);
				// continue;
			}

			else
			{
				p3.push_back(chain2[j][0] + offx);
				p3.push_back(chain2[j][1] + offy);
				p3.push_back(chain2[j][2] + offz);
				p4.push_back(chain2[j + 1][0] + offx);
				p4.push_back(chain2[j + 1][1] + offy);
				p4.push_back(chain2[j + 1][2] + offz);
			}

			result += compute_one(p1, p2, p3, p4);
		}
	}

	return result;
}

/**
 * Compute the inverse of a given matrix
 *
 * @param m The matrix to invert
 * @return The inverse of the matrix m
 */
vector<vector<double>> invert_mat(vector<vector<double>> m)
{
	vector<vector<double>> ans;
	vector<double> temp1;
	vector<double> temp2;
	vector<double> temp3;
	double den = m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1];
	den += -m[1][0] * m[0][1] * m[2][2] + m[1][0] * m[0][2] * m[2][1];
	den += m[2][0] * m[0][1] * m[1][2] - m[2][0] * m[0][2] * m[1][1];

	temp1.push_back((m[1][1] * m[2][2] - m[1][2] * m[2][1]) / den);
	temp1.push_back(-(m[0][1] * m[2][2] - m[0][2] * m[2][1]) / den);
	temp1.push_back((m[0][1] * m[1][2] - m[0][2] * m[1][1]) / den);
	temp2.push_back(-(m[1][0] * m[2][2] - m[1][2] * m[2][0]) / den);
	temp2.push_back((m[0][0] * m[2][2] - m[0][2] * m[2][0]) / den);
	temp2.push_back(-(m[0][0] * m[1][2] - m[0][2] * m[1][0]) / den);
	temp3.push_back((m[1][0] * m[2][1] - m[1][1] * m[2][0]) / den);
	temp3.push_back(-(m[0][0] * m[2][1] - m[0][1] * m[2][0]) / den);
	temp3.push_back((m[0][0] * m[1][1] - m[0][1] * m[1][0]) / den);

	ans.push_back(temp1);
	ans.push_back(temp2);
	ans.push_back(temp3);

	return ans;
}

/**
 * Generate a matrix that can be multiplied by coordinates to randomly project those coordinates
 *
 * @return A matrix to be used to generate random projections
 */
vector<vector<double>> get_random_proj()
{
	int error;
	ruv_buffer_t buffer;
	int ret;
	vector<vector<double>> matrix;
	double* rvector = new double[3];

	ret = ruv_init(&buffer);
	ret = ruv_generate(&buffer, 3, rvector);
	vector<double> temp;
	for (int i = 0; i < 3; i++)
	{
		temp.push_back(rvector[i]);
	}
	
	matrix.push_back(temp);
	while (1)
	{
		double vec[3];
		double* cross = new double[3];
		ruv_generate(&buffer, 3, vec);
		cross3(rvector, vec, cross);
		if (abs(cross[0]) > 0.0001 || abs(cross[1]) > 0.0001 || abs(cross[2]) > 0.0001)
		{
			
			vector<double> temp1;
			for (int i = 0; i < 3; i++)
			{
				temp1.push_back(cross[i]);
			}
			matrix.push_back(temp1);
			delete[] cross;
			break;
		}
		else
		{
			delete[] cross;
		}
	}

	while (1)
	{
		double vec[3];
		double* cross = new double[3];
		ruv_generate(&buffer, 3, vec);
		cross3(rvector, vec, cross);
		if (abs(cross[0]) > 0.0001 || abs(cross[1]) > 0.0001 || abs(cross[2]) > 0.0001)
		{
			
			vector<double> temp2;
			for (int i = 0; i < 3; i++)
			{
				temp2.push_back(cross[i]);
			}
			matrix.push_back(temp2);
			delete[] cross;
			break;
		}
		else
		{
			delete[] cross;
		}
	}

	ruv_free(&buffer);
	delete[] rvector;
	double trans[3][3];
	vector<vector<double>> transposed;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			trans[i][j] = matrix[j][i];
		}
		transposed.push_back({trans[i][0], trans[i][1], trans[i][2]});
	}
	vector<vector<double>> inverted = invert_mat(transposed);

	return inverted;
}

/**
 * Unwrap the coordinates of a point within a system that uses periodic boundary conditions
 *
 * @param p1, p2 Two adjacent points in a chain
 * @param box_dim The dimensions of the periodic box for the given periodic system
 * @return The unwrapped coordinates for p2
 */
vector<double> unwrap(vector<double> p1, vector<double> p2, double box_dim)
{
	int imgx = round((p1[0] - p2[0]) / box_dim);
	int imgy = round((p1[1] - p2[1]) / box_dim);
	int imgz = round((p1[2] - p2[2]) / box_dim);

	return {p2[0] + (imgx * box_dim), p2[1] + (imgy * box_dim), p2[2] + (imgz * box_dim)};
}

/**
 * Obtain the point of intersection between two edges
 *
 * @param p0, p1 The points that form the first edge in the calculation
 * @param p2, p3 The points that form the second edge in the calculation
 * @return The point of intersection between the two edges
 */
vector<double> get_intersection(vector<double> p0, vector<double> p1, vector<double> p2, vector<double> p3)
{
	double s1_x, s1_y, s2_x, s2_y, det, t;
	s1_x = p1[0] - p0[0];
	s1_y = p1[1] - p0[1];
	s2_x = p3[0] - p2[0];
	s2_y = p3[1] - p2[1];
	det = -s2_x * s1_y + s1_x * s2_y;
	t = (s2_x * (p0[1] - p2[1]) - s2_y * (p0[0] - p2[0])) / det;
	vector<double> result = {p0[0] + (t * s1_x), p0[1] + (t * s1_y)};
	return result;
}

/**
 * Multiply two polynomials, used in calculating the Jones polynomial
 *
 * @param power The number of separate loops in the system being measured
 * @return The polynomial of a single combination of crossing annealments
 */
map<int, double> mult_poly(int power)
{
	map<int, double> base;
	base[-2] = -1;
	base[2] = -1;
	if (power == 2)
	{
		return base;
	}

	map<int, double> result;
	map<int, double> temp;
	temp[-2] = -1;
	temp[2] = -1;
	for (int k = 2; k < power; k++)
	{
		for (auto const& x: temp)
		{
			for (auto const& y: base)
			{
				result[x.first + y.first] += x.second * y.second;
			}
		}
		temp = result;
		result.clear();
	}

	return temp;
}

/**
 * Multiply two polynomials where both are passed as parameters
 *
 * @param a, b The two polynomials to multiply
 * @return The product of a and b
 */
map<int, double> simple_mult(map<int, double> a, map<int, double> b)
{
	map<int, double> result;
	for (auto const& x: a)
	{
		for (auto const& y: b)
		{
			result[x.first + y.first] += x.second * y.second;
		}
	}

	return result;
}

/**
 * Create a vector of vectors that each contain the indices of a given atom's neighbor atoms.
 * Ex: neigh_array[1] contains the neighbors of atom 1, which are 0 and 2.
 *
 * @param n The number of atoms in the system
 * @param closed True if the last atom in the system is connected to the first, false otherwise
 * @return A vector of vectors containing each atom's neighbors
 */
vector<vector<int>> generate_neigh_array(int n, bool closed)
{
	vector<vector<int>> neigh_array;
	for (int i = 1; i < n - 1; i++)
	{
		neigh_array.push_back({i - 1, i + 1});
	}

	if (closed)
	{
		neigh_array.push_back({n - 2, 0});
		neigh_array.insert(neigh_array.begin(), {n - 1, 1});
	}
	else
	{
		neigh_array.push_back({n - 2, -1});
		neigh_array.insert(neigh_array.begin(), {-1, 1});
	}

	return neigh_array;
}

/**
 * Finds all the cells that the parent image of a chain intersects.
 * Cells are represented by a vector containing the x, y, and z image flags.
 * The final results contains a vector of all the cells that the parent image intersects.
 *
 * @param coords The coordinates of the system to find the images of
 * @param n The number of atoms in the system
 * @param box_dims The dimensions of the periodic box
 * @return A vector of vectors each containing a unique image used by the given system
 */
vector<vector<int>> compute_img(double** coords, int n, vector<double> box_dims)
{
	vector<vector<int>> result;
	for (int i = 0; i < n; i++)
	{
		/**
		 * We assume input chains (unfolded) in an orgin centered box
		 *
		 */
		int ximg = round(coords[i][0] / (box_dims[0]));
		int yimg = round(coords[i][1] / (box_dims[1]));
		int zimg = round(coords[i][2] / (box_dims[2]));
		vector<int> temp = {ximg, yimg, zimg};
		bool found = false;
		for (auto vec: result)
		{
			if (vec == temp)
			{
				found = true;
				break;
			}
		}
		if (!found)
		{
			result.push_back(temp);
		}
	}

	return result;
}




/**
 * Finds all the images that intersect all the cells that the parent image of a chain intersects.
 * Images are represented by a vector containing the x, y, and z image flags.
 * The final results contains a vector of all the images that intersect the cells in which the parent image lies in.
 *
 * @param cells The cells that the parent image of a chain intersects
 * @param coords The coordinates of the system to find the images of
 * @param n The number of atoms in the system
 * @param box_dims The dimensions of the periodic box
 * @return A vector of vectors each containing a unique image used by the given system
 */



vector<vector<int>> compute_periodic_img(vector<vector<int>> cells, double** coords, int n, vector<double> box_dims)
{
        vector<vector<int>> result;
        for (int i = 0; i < cells.size()-1; i++)
        {
            for(int j=i+1; j<cells.size();j++)   
                {
                /**
                 * We assume input chains (unfolded) in an orgin centered box
                 *
                 */
                int ximg = -1*cells[i][0]+cells[j][0];
                int yimg = -1*cells[i][1]+cells[j][1];
                int zimg = -1*cells[i][2]+cells[j][2];
                vector<int> temp = {ximg, yimg, zimg};
                bool found = false;
                for (auto vec: result)
                {
                        if (vec == temp)
                        {
                                found = true;
                                break;
                        }
                }
                if (!found)
                {
                        result.push_back(temp);
                }
                }
        }

        return result;
}




/**
 * Translates a chain by a given vector
 * The final result is the coordinates of the translated chain.
 *
 * @param coords The coordinates of a chain
 * @param n The number of atoms in the system
 * @param translation_vector The vector by which we translate the chain
 * @param box_dims The dimensions of the periodic box
 * @return The coordinates of the translated chain
 */



//vector<vector<double>> 
double** translation(double** coords,int n,  vector<int> translation_vector , vector<double> box_dims)
{
  //vector<vector<double>> result;
  vector<vector<double>> temp_result;
  for (int i=0; i<n;i++)
  {
     double ximg = (coords[i][0]+translation_vector[0])*box_dims[0];
     double yimg = (coords[i][1]+translation_vector[1])*box_dims[1];
     double zimg = (coords[i][2]+translation_vector[2])*box_dims[2];
     vector<double> temp = {ximg, yimg, zimg};
     temp_result.push_back(temp);
  }
  double** result = new double*[temp_result.size()];
  for (int i = 0; i < temp_result.size(); i++)
  {
       result[i] = new double[3];
       copy(temp_result[i].begin(), temp_result[i].end(), result[i]);
  }
   
  return result;

}

/**
 * Counts the number of crossings in a projection
 *
 * @param neigh_array A vector of vectors containing each atom's neighboring atoms
 * @param proj The projected coordinates of the system to count the crossings of
 * @param n The number of atoms in the system
 * @param count_all True if the last atom in the system is connected to the first, false otherwise
 * @return A vector of vectors each containing 4 atoms that form two intersecting edges
 */
vector<vector<int>> count_crossings(vector<vector<int>> neigh_array, vector<vector<double>> proj, int n, bool is_closed)
{
	int count = 0;
	vector<vector<int>> res;
        vector<vector<int>> final_res;
	int lasti, lastj;
        lasti = n-2;
	map<int, vector<vector<int>>> order;
	map<int, vector<double>> distances;

	for (int i = 0; i < lasti; i++)
	{
		if (neigh_array[i][1] < 0)
			continue;
                lastj = n;
		for (int j = i + 2; j < lastj; j++)
		{
			int p2 = j % n;
			if (abs(i - p2) < 2)
				continue;
			if (i == neigh_array[p2][1])
				continue;
			if (neigh_array[i][1] == p2)
				continue;
			if (neigh_array[p2][1] < 0)
				continue;
			double rx = 1;
			double ry = 1;
			if (intersect1(proj[i], proj[neigh_array[i][1]], proj[p2], proj[neigh_array[p2][1]], &rx, &ry))
			{
				double dist = distance(proj[i], {rx, ry});
				if (order.find(i) == order.end())
				{
					order[i].push_back({i, neigh_array[i][1], p2, neigh_array[p2][1]});
					distances[i].push_back(dist);
				}
				else
				{
					bool found = false;
					for (int z = 0; z < order[i].size(); z++)
					{
						if (dist < distances[i][z])
						{
							order[i].insert(order[i].begin() + z, {i, neigh_array[i][1], p2, neigh_array[p2][1]});
							distances[i].insert(distances[i].begin() + z, dist);
							found = true;
							break;
						}
					}
					if (!found)
					{
						order[i].push_back({i, neigh_array[i][1], p2, neigh_array[p2][1]});
						distances[i].push_back(dist);
					}
				}
				if (i == neigh_array[p2][1])
					continue;
				if (neigh_array[i][1] == p2)
					continue;
			}
		}
	}
	for (auto it = order.begin(); it != order.end(); it++)
	{
		int i = it->first;
		for (int j = 0; j < order[i].size(); j++)
		{
			res.push_back(order[i][j]);
		}
	}
	return res;
}

/**
 * Reduces the number of crossings in a system by performing Reidemeister I and II moves.
 *
 * @param initial_crossings A vector of vectors each containing 4 atoms that form two intersecting edges
 * @param coords The coordinates of the system
 * @return A reduced vector of vectors each containing 4 atoms that form two intersecting edges
 */
vector<vector<int>> reduce_crossings(vector<vector<int>> initial_crossings, vector<vector<double>> coords)
{
	int init_crossings = initial_crossings.size();
	map<int, int> crossing_map;
	for (int i = 0; i < init_crossings; i++)
	{
		crossing_map[initial_crossings[i][0]] = initial_crossings[i][1];
		crossing_map[initial_crossings[i][2]] = initial_crossings[i][3];
	}
	vector<vector<int>> temp_result;
	vector<int> remove_crossings;
	int count = 0;
	bool boolarray[init_crossings];  
	for (int i = 0; i < init_crossings; i++)
	{
		boolarray[i] = true;
	}
	int p11, p12, p13, p14, p21, p22, p23, p24, p31, p32, p33, p34;
	for (int i = 0; i < initial_crossings.size(); i++)
	{
                
		bool found = false;
		if (find(remove_crossings.begin(), remove_crossings.end(), i) != remove_crossings.end())
			continue;
		p11 = initial_crossings[i][0];
		p12 = initial_crossings[i][1];
		p13 = initial_crossings[i][2];
		p14 = initial_crossings[i][3];
		if (found)
			continue;
		double wr1 = compute_one(
			coords[initial_crossings[i][0]],
			coords[initial_crossings[i][1]],
			coords[initial_crossings[i][2]],
			coords[initial_crossings[i][3]]);
		for (int j = i + 1; j < initial_crossings.size(); j++)
		{
                        
			found = false;
			if (find(remove_crossings.begin(), remove_crossings.end(), j) != remove_crossings.end())
				continue;
			if (found)
				continue;
			double wr2 = compute_one(
				coords[initial_crossings[j][0]],
				coords[initial_crossings[j][1]],
				coords[initial_crossings[j][2]],
				coords[initial_crossings[j][3]]);
                        
			if ((wr1 < 0 && wr2 < 0) || wr1 > 0 && wr2 > 0)
				found = true;
                                
			if (found)
				continue;
			p21 = initial_crossings[j][0];
			p22 = initial_crossings[j][1];
			p23 = initial_crossings[j][2];
			p24 = initial_crossings[j][3];
          
			if (p12 == p21 && (p14 == p23 || p24 == p13))
			{
				remove_crossings.push_back(i);
				remove_crossings.push_back(j);
				j = initial_crossings.size();
				continue;
			}
			if (found)
				continue;
		}
	}

	for (int i = remove_crossings.size() - 1; i >= 0; i--)
	{
		initial_crossings.erase(initial_crossings.begin() + remove_crossings[i]);
	}
	remove_crossings.clear();
	int p1, p2, p3, p4;
	for (int i = 0; i < initial_crossings.size(); i++)
	{
		bool found = false;
		p1 = initial_crossings[i][0];
		p2 = initial_crossings[i][1];
		p3 = initial_crossings[i][2];
		p4 = initial_crossings[i][3];
		if (p2 == 0)
			p2 += coords.size();
		if (p4 == 0)
			p4 += coords.size();
		for (int j = min(p2, p3); j < max(p2, p3); j++)
		{
			if (crossing_map.find(j) != crossing_map.end())
			{
				found = true;
				break;
			}
			if (found)
				continue;
		}
		if (!found)
		{
			remove_crossings.push_back(i);
			crossing_map.erase(initial_crossings[i][0]);
			crossing_map.erase(initial_crossings[i][2]);
		}
	}

	sort(remove_crossings.begin(), remove_crossings.end());
	for (int i = remove_crossings.size() - 1; i >= 0; i--)
	{
		initial_crossings.erase(initial_crossings.begin() + remove_crossings[i]);
	}

	count = 0;
	vector<vector<int>> result;
	for (int i = 0; i < initial_crossings.size(); i++)
	{
		result.push_back(
			{initial_crossings[i][0], initial_crossings[i][1], initial_crossings[i][2], initial_crossings[i][3]});
	}

	return result;
}

/**
 * Reads unwrapped coordinates from a .tepp file.
 *
 * @param[in] fname The name (including path) of the .tepp file
 * @param[out] n The number of atoms in the .tepp file
 * @return A 2d array containing the coordinates for each atom in the system
 */
double** read_coords(string fname, int* n)
{
	string line;
	ifstream myfile;
	myfile.open(fname);
	int num_atoms = 0;
	string::size_type sz;
	vector<vector<double>> temp;

	if (!myfile.is_open())
	{
		cout << "Couldn't open file\n";
		exit(EXIT_FAILURE);
	}

	while (getline(myfile, line))
	{
		string buf;
		stringstream ss(line);
		int count = 0;
		vector<string> tokens;
		while (ss >> buf)
		{
			tokens.push_back(buf);
			count++;
		}

		if (count > 0)
		{
			temp.push_back({stod(tokens[0], &sz), stod(tokens[1], &sz), stod(tokens[2], &sz)});
		}
	}

	double** result = new double*[temp.size()];
	for (int i = 0; i < temp.size(); i++)
	{
		result[i] = new double[3];
		copy(temp[i].begin(), temp[i].end(), result[i]);
	}

	*n = temp.size();
	return result;
}

/**
 * Reads wrapped coordinates from a .tepp file and unwraps them.
 *
 * @param[in] fname The name (including path) of the .tepp file
 * @param[out] n The number of atoms in the .tepp file
 * @param[in] chain_length The length of each individual chain within the system
 * @param[in] box_dim The dimensions of the periodic box
 * @return A 2d array containing the unwrapped coordinates for each atom in the system
 */
double** read_coords(string fname, int* n, int chain_length, double box_dim)
{
	string line;
	ifstream myfile;
	myfile.open(fname);
	int num_atoms = 0;
	string::size_type sz;
	vector<vector<double>> temp;

	if (!myfile.is_open())
	{
		cout << "Couldn't open file\n";
		exit(EXIT_FAILURE);
	}

	while (getline(myfile, line))
	{
		string buf;
		stringstream ss(line);
		int count = 0;
		vector<string> tokens;
		while (ss >> buf)
		{
			tokens.push_back(buf);
			count++;
		}

		if (count > 0)
		{
			if (temp.size() == 0 || (temp.size() % chain_length) == 0)
			{
				temp.push_back({stod(tokens[0], &sz), stod(tokens[1], &sz), stod(tokens[2], &sz)});
			}
			else
			{
				temp.push_back(unwrap(
					temp[temp.size() - 1],
					{stod(tokens[0], &sz), stod(tokens[1], &sz), stod(tokens[2], &sz)},
					box_dim));
			}
		}
	}

	double** result = new double*[temp.size()];
	for (int i = 0; i < temp.size(); i++)
	{
		result[i] = new double[3];
		copy(temp[i].begin(), temp[i].end(), result[i]);
	}

	*n = temp.size();
	return result;
}

/**
 * Projects the coordinates of a system using random projection vectors.
 *
 * @param coords The coordinates of the system to project
 * @param n The number of atoms in the system
 * @return The projected coordinates of the system
 */
vector<vector<double>> get_proj(vector<vector<double>> coords, int n)
{
	vector<vector<double>> result;
	vector<vector<double>> matrix = get_random_proj();

	for (int i = 0; i < n; i++)
	{
		double xk
			= ((matrix[0][0] * coords[i][0]) + (matrix[0][1] * coords[i][1]) + (matrix[0][2] * coords[i][2]))
			  / sqrt((matrix[0][0] * matrix[0][0]) + (matrix[0][1] * matrix[0][1] + (matrix[0][2] * matrix[0][2])));
		double yk
			= ((matrix[1][0] * coords[i][0]) + (matrix[1][1] * coords[i][1]) + (matrix[1][2] * coords[i][2]))
			  / sqrt((matrix[1][0] * matrix[1][0]) + (matrix[1][1] * matrix[1][1] + (matrix[1][2] * matrix[1][2])));
		double zk
			= ((matrix[2][0] * coords[i][0]) + (matrix[2][1] * coords[i][1]) + (matrix[2][2] * coords[i][2]))
			  / sqrt((matrix[2][0] * matrix[2][0]) + (matrix[2][1] * matrix[2][1] + (matrix[2][2] * matrix[2][2])));
		vector<double> temp = {xk, yk, zk};
		result.push_back(temp);
	}

	return result;
}
 

/**
 * Finds edges that have intersect multiple other edges and breaks such edges so that
 * every edge has at most one intersection.
 *
 * @param proj The projected coordinates of the system
 * @param coords The actual coordinates of the system
 * @param neigh_array A vector of vectors containing each atom's neighboring atoms
 * @param n The number of atoms in the system
 * @param is_closed True if the last atom in the system is connected to the first, false otherwise
 * @return A struct containing the new projected coordinates, actual coordinates, and vector of crossings after all
 * edges with multiple crossings have been split
 */
Struct mult_crossings(
        vector<vector<double>> proj, vector<vector<double>> coords, vector<vector<int>> neigh_array, vector<vector<int>> before_cross, int n, bool is_closed)
{
        int last = -1;
        int off = 0;
        int last_off = 0;
        int temp_off;
        map<int, int> offmap;
        map<int, int> same_crossing;
        int offset=0;
        vector<vector<double>> new_proj;
        vector<vector<double>> new_coords; 
        vector<vector<int>> new_before_cross;
        int lastedge=n;
        
        
        map<int, vector<int>> multiple_crossings = has_mult_crossings(n,before_cross,is_closed); 
         
        vector<vector<int>> new_cross; 
        vector<int> splitat;    
        vector<int> countm;
        for (int k=0;k<lastedge;k++)
        {
           
            vector<double> fcr;
            vector<double> distances;
            new_proj.push_back(proj[k]);
            new_coords.push_back(coords[k]);
            if (multiple_crossings[k].size()<2)
               countm.push_back(1);
            if (multiple_crossings[k].size()>1)
            {
               
               int ls=multiple_crossings.size();
               countm.push_back(ls);
               vector<vector<double>> crossings;
               int next_p = k + 1;
               if (k==n-1)
                   next_p=0;
               for (int s=0; s < multiple_crossings[k].size();s++)
               {
                   
                   int lv=multiple_crossings[k][s];
                   int lvn;
                   if (lv==n-1)
                   {
                      lvn=0;
                   }
                   else
                   {
                      lvn=lv+1;
                   }
                   
                   vector<double> temp = get_intersection( proj[k], proj[next_p], proj[lv], proj[lvn]);
                   crossings.push_back(temp);
                   double dist=distance(proj[k],temp);
                   
                   if (s==0)
                   {
                      
                      fcr.insert(fcr.begin()+0,0);
                      distances.insert(distances.begin()+0,dist);
                      
                   }
                   else
                   {
                      if (dist<distances[0])
                      {
                         auto a=0;
                         
                         fcr.insert(fcr.begin()+0,s);
                         distances.insert(distances.begin()+0,dist);
                      }
                      else
                      {
                         if (dist>distances[distances.size()-1])
                         {
                            int r=distances.size();
                            
                            fcr.insert(fcr.begin()+r, s);
                            distances.insert(distances.begin()+r,dist);
                         }
                         else
                         {
                            for (int h=1;h<distances.size();h++)
                            {
                           
                               auto a=h;
                               if ((dist<=distances[h]) and (dist>distances[h-1]))
                               {
                                  
                                  fcr.insert(fcr.begin()+h, s);
                                  distances.insert(distances.begin()+h,dist);
                                  
                                  break;
                               }
                            }
                         }
                      }
                   }
                   
                   
                   
               }
               
               
               for (int t=1; t < multiple_crossings[k].size();t++)
               {       
                   int s = fcr[t];
                   int u = fcr[t-1];
                   
                   double newx2 = crossings[u][0] + ((crossings[s][0] - crossings[u][0]) / 2);
                   double newy2 = crossings[u][1] + ((crossings[s][1] - crossings[u][1]) / 2);
                   double newx_coords
                                        = coords[k][0] + (t) * ((coords[next_p][0] - coords[k][0]) / multiple_crossings[k].size());
                   double newy_coords
                                        = coords[k][1] + (t) * ((coords[next_p][1] - coords[k][1]) / multiple_crossings[k].size());
                   double newz_coords
                                        = coords[k][2] + (t) * ((coords[next_p][2] - coords[k][2]) / multiple_crossings[k].size());
                   
                      
                   vector<double> pnt = {newx2, newy2};
                   new_proj.push_back(pnt);
                   new_coords.push_back({newx_coords, newy_coords, newz_coords});
                   
                               
               }
            }
        }
        
        
        Struct final;
        final.coords = new_coords;
        final.proj = new_proj;
        final.success = true;
        return final;

}





/**
 * Checks the given polynomial to see if it matches the polynomial of the given knot type(s).
 *
 * @param poly The polynomial to check
 * @param params A vector containing the names of knots to check poly against for matches
 * @return A vector with the same length as params with a 1 at result[i] if poly matches the knot type at params[i], or
 * a 0 at result[i] otherwise
 */
vector<int> find_knot(map<int, double> poly, vector<string> params)
{
	vector<int> result;
	map<string, map<int, double>> knots;
	map<int, double> trefoil1 = {{-16, -1}, {-12, 1}, {-4, 1}};
	knots["trefoil1"] = trefoil1;
	map<int, double> trefoil2 = {{4, 1}, {12, 1}, {16, -1}};
	knots["trefoil2"] = trefoil2;
	map<int, double> figure8 = {{-8, 1}, {-4, -1}, {0, 1}, {4, -1}, {8, 1}};
	knots["figure8"] = figure8;
	map<int, double> pentafoil = {{-28, -1}, {-24, 1}, {-20, -1}, {-16, 1}, {-8, 1}};
	knots["pentafoil"] = pentafoil;
	map<int, double> stevedore = {{-16, 1}, {-12, -1}, {-8, 1}, {-4, -2}, {0, 2}, {4, -1}, {8, 1}};
	knots["stevedore"] = stevedore;
	for (int i = 0; i < params.size(); i++)
	{
		if (params[i].compare("trefoil") == 0)
		{
			bool found = false;
			for (int j = 0; j < knots["trefoil1"].size(); j++)
			{
				if (poly.find(knots["trefoil1"][j]) == poly.end() || abs(poly[knots["trefoil1"][j]]) < 0.1)
				{
					found = true;
					break;
				}
			}
			if (!found)
				result.push_back(1);
			found = false;
			for (int j = 0; j < knots["trefoil2"].size(); j++)
			{
				if (poly.find(knots["trefoil2"][j]) == poly.end() || abs(poly[knots["trefoil1"][j]]) < 0.1)
				{
					found = true;
					break;
				}
			}
			if (!found)
				result.push_back(1);
			else
				result.push_back(0);
		}
		else if (params[i].compare("figure8") == 0)
		{
			bool found = false;
			for (int j = 0; j < knots["figure8"].size(); j++)
			{
				if (poly.find(knots["figure8"][j]) == poly.end() || abs(poly[knots["figure8"][j]]) < 0.1)
				{
					found = true;
					break;
				}
			}
			if (!found)
				result.push_back(1);
			else
				result.push_back(0);
		}
		else if (params[i].compare("pentafoil") == 0)
		{
			bool found = false;
			for (int j = 0; j < knots["pentafoil"].size(); j++)
			{
				if (poly.find(knots["pentafoil"][j]) == poly.end() || abs(poly[knots["pentafoil"][j]]) < 0.1)
				{
					found = true;
					break;
				}
			}
			if (!found)
				result.push_back(1);
			else
				result.push_back(0);
		}
		else if (params[i].compare("stevedore") == 0)
		{
			bool found = false;
			for (int j = 0; j < knots["stevedore"].size(); j++)
			{
				if (poly.find(knots["stevedore"][j]) == poly.end() || abs(poly[knots["stevedore"][j]]) < 0.1)
				{
					found = true;
					break;
				}
			}
			if (!found)
				result.push_back(1);
			else
				result.push_back(0);
		}
		else
		{
			cout << params[i] << " knot not supported!\n";
		}
	}

	return result;
}

/**
 * This function creates an output directory if it doesn't exist,
 * if it does exists, skip creating the directory.
 *
 */
void create_output_dir()
{
	if (!fs::exists("output"))
	{
		cout << "Creating output directory..." << endl;
		fs::create_directory("output");
	}
}

#endif
