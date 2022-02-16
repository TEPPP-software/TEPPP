#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
namespace fs = std::filesystem;

/**
 * Converts a file with extension .pdb, .dcd, or .read_data
 * into a .teppp file that can be read by the rest of the functions.
 * Converted files are placed in $INSTALL_DIR/converted.
 *
 * @param[in] argv[1] Name (including path) of file to be converted
 */
int main(int argc, char** argv)
{
	// Check if converted directory exists, if not create it
	if (!fs::exists("converted"))
	{
		std::cout << "Creating converted directory" << std::endl;
		fs::create_directory("converted");
	}
	else
	{
		std::cout << "converted directory already exists. Skipping..." << std::endl;
	}
	std::string infname(argv[1]);
	std::string outfname = "./converted/" + std::string(fs::path(infname).filename());
	std::string ext = std::string(fs::path(infname).extension());
	fs::path p = fs::path(outfname).replace_extension(".teppp");
	outfname = std::string(p);
	std::ofstream outfile;
	outfile.open(outfname);

	if (ext.compare(".read_data") == 0)
	{
		std::string line;
		std::ifstream myfile;
		myfile.open(infname);
		int num_atoms = 0;
		std::string::size_type sz;

		if (!myfile.is_open())
		{
			std::cout << "Couldn't open file\n";
			exit(EXIT_FAILURE);
		}

		while (getline(myfile, line))
		{
			std::string buf;
			std::stringstream ss(line);
			int count = 0;
			std::vector<std::string> tokens;
			while (ss >> buf)
			{
				tokens.push_back(buf);
				count++;
			}

			if (count == 2 && tokens[1].compare("atoms") == 0)
			{
				num_atoms = stoi(tokens[0], &sz);
				break;
			}
		}

		if (num_atoms == 0)
		{
			std::cout << "num_atoms not found\n";
			exit(EXIT_FAILURE);
		}

		while (getline(myfile, line))
		{
			if (line.compare("Atoms") == 0)
			{
				getline(myfile, line);
				for (int i = 0; i < num_atoms; i++)
				{
					getline(myfile, line);
					std::string buf;
					std::stringstream ss(line);
					int count = 0;
					std::vector<std::string> tokens;
					while (ss >> buf)
					{
						tokens.push_back(buf);
						count++;
					}

					if (count > 0)
					{
						outfile << tokens[3] + " " + tokens[4] + " " + tokens[5] + "\n";
					}
				}
				break;
			}
		}

		myfile.close();
		outfile.close();
	}
	else
	{
		std::cout << "File extension not supported. Exiting...\n";
		exit(EXIT_FAILURE);
	}
	return 0;
}
