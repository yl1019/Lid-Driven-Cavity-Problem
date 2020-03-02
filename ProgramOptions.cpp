#include<iostream>
#include<cstdlib>
using namespace std;

#include <boost/program_options.hpp>

namespace po = boost::program_options;

/**
 * @brief Validate if the input from command line is read successfully
 * @param argc
 * @param argv
 * @param vm a map to store the command-line arguments
 * @return Return 0 if help is called or an error occured; return 1 if parameters from command line are read successfully and the main program continue to proceed.
 */
bool OptionStatus(int argc, char *argv[], po::variables_map &vm)
{
   try
   {
	po::options_description desc("Allowed input options");
	/// Adding input options with appropriate default values
	desc.add_options()
		("help", "Produce help message")
		("Lx", po::value<double>()->default_value(1.0), "Length of the domain in the x-direction")
		("Ly", po::value<double>()->default_value(1.0), "Length of the domain in the y-direction")
		("Nx", po::value<int>() ->default_value(161), "Number of grid points in the x-direction")
		("Ny", po::value<int>() ->default_value(161), "Number of grid points in the y-direction")
		("Px", po::value<int>() ->default_value(1), "Number of the partitions in the x-direction")
		("Py", po::value<int>() ->default_value(1), "Number of the partitions in the y-direction")
		("dt", po::value<double>() ->default_value(0.01),"Time step size")
		("T", po::value<double>() ->default_value(5.0),"Final time")
		("Re", po::value<double>() ->default_value(100.0),"Reynolds")
	;
	/// Parse the command-line arguments and store in buffer vm
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	/// If user asked for help, display help information and return 0
	if (vm.count("help"))
	{
		cout << "Usage: option_description [options]\n";
		cout << desc << endl;
		return false;
	}
   }

   catch(exception const &e)
   {
	cout << e.what() << endl;
	return false;
   }
   return true;
}
/**
 * @brief collect datas from the buffer vm
 */
void ReadVals(po::variables_map &vm, double &Lx, double &Ly, int &Nx, int &Ny, int &Px, int &Py, double &dt, double &T, double &Re)
{
	Lx = vm["Lx"].as<double>();
	Ly = vm["Ly"].as<double>();
	Nx = vm["Nx"].as<int>();
	Ny = vm["Ny"].as<int>();
	Px = vm["Px"].as<int>();
	Py = vm["Py"].as<int>();
	dt = vm["dt"].as<double>();
	T = vm["T"].as<double>();
	Re = vm["Re"].as<double>();
}
