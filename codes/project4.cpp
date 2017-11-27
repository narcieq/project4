// project4.cpp : 콘솔 응용 프로그램에 대한 진입점을 정의합니다.
//
/*
#include "mpi.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>
using namespace std;

ofstream ofile;


int main(int argc, char* argv[])
{
	/*MPI_Init(&argc, &argv);//initailize mpi status by command line

	int rank, numprocs;

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);// store total number of pc into "numprocs"
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // store current pc rank into "rank"

	if (rank == 0 && argc <= 1) {
		cout << "Bad Usage: " << argv[0] << " read output file" << endl;
		exit(1);
	}
	else if (rank == 0 && argc > 1) {
		string outfilename;
		outfilename = argv[1];
		ofile.open(outfilename);
	}
	cout << numprocs << " " << rank << endl;

	MPI_Finalize();

	//set lattice
	int L = 2;
	double** lattice;
	lattice = new double*[L];
	for (int i = 0; i < L; i++) {
		lattice[L] = new double[L];
		for (int j = 0; j < L; j++) {
			lattice[L][j] = 0;
		}
	}

	srand((unsigned int)time(0));

	int random = 0;
	double M = 0;

	//initialize random spin to lattice
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			do {
				random = rand() % (2 + 1) - 1;
			} while (random == 0);
			lattice[i][j] = random;
			M += lattice[i][j];
		}
	}

	//calculate Ei;
	double* E = new double[L*L];
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			E[i + j] = lattice[i][j] * (lattice[p][j] + lattice[i][p]);
		}
	}

	//filp spin i
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			lattice[i][j] = -lattice[i][j];
		}
	}

	return 0;
}

int boundary(int index, int lattice_size, int add) {
	return (index + lattice_size + add) % lattice_size;
}
*/



#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <vector>
//#include <armadillo>
#include <string>
using namespace  std;
//using namespace arma;
// output file
ofstream ofile;

double ran2(long *);

//structure to count probability
struct e_count {
	double eng;
	int count = 0;
};
vector<e_count> E_vector;
vector<e_count>E_total_vector;

// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
	return (i + limit + add) % (limit);
}
// Function to initialise energy and magnetization
void InitializeLattice(int, double**, double&, double&);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, double*, long&, int, int);
// prints to file the results of the calculations  
void WriteResultstoFile(int, int, double, double*);

// Main program begins here

int main(int argc, char* argv[])
{
	string filename;
	int NSpins, MCcycles;
	double InitialTemp, FinalTemp, TempStep;
	if (argc <= 5) {
		cout << "Bad Usage: " << argv[0] <<
			" read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
		exit(1);
	}
	if (argc > 1) {
		filename = argv[1];
		NSpins = atoi(argv[2]);
		MCcycles = atoi(argv[3]);
		InitialTemp = atof(argv[4]);
		FinalTemp = atof(argv[5]);
		TempStep = atof(argv[6]);
	}
	// Declare new file name and add lattice size to file name
	string fileout = filename;
	string argument = to_string(NSpins);
	fileout.append(argument);
	argument = to_string(MCcycles);
	fileout += "_";
	fileout.append(argument);
	fileout += ".txt";
	ofile.open(fileout);

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

	// Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);



	cout << "processor name : " << processor_name << endl;
	cout << "rank : " << world_rank << endl;
	cout << "size : " << world_size << endl;

	//#pragma mpi parallel for default(shared) private(NSpins, MCcycles, Temperature, ExpectationValues) reduction()
		// Start Monte Carlo sampling by looping over the selcted Temperatures
	for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature += TempStep) {
		double* ExpectationValues = new double[6];
		double TotalExpectation[6];
		for (int i = 0; i < 6; i++) {
			ExpectationValues[i] = 0;
			TotalExpectation[i] = 0;
		}

		int no_intervalls = MCcycles / world_size;
		int myloop_begin = world_rank*no_intervalls + 1;
		int myloop_end = (world_rank + 1)*no_intervalls;
		if ((world_rank == world_size - 1) && (myloop_end < MCcycles)) myloop_end = MCcycles;

		// Start Monte Carlo computation and get expectation values
		MPI_Bcast(&NSpins, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&InitialTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&FinalTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&TempStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		long idum = -1 - world_rank;

		double  TimeStart, TimeEnd, TotalTime;
		TimeStart = MPI_Wtime();

		MetropolisSampling(NSpins, MCcycles, Temperature, ExpectationValues, idum, myloop_begin, myloop_end);

		cout << "MCS finished" << endl;
		
		for (int i = 0; i < 6; i++) {
			int a = MPI_Reduce(&ExpectationValues[i], &TotalExpectation[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			
		}
/*
		for (int i = 0; i < (int)E_vector.size(); i++) {
			e_count temp;
			temp.count = 0;
			temp.eng = 0;
			E_total_vector.push_back(temp);
			for (int j = 0; j < (int)E_vector.size(); j++) {
				int* sameE = new int;
				if (MPI_Comm_compare(0, 1, sameE) == MPI_IDENT) {
					MPI_Reduce(&E_vector[i].count, &E_total_vector[i].count, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				}
			}
		}
*/		
		// print results
		if (world_rank == 0) {
			WriteResultstoFile(NSpins, MCcycles, Temperature, TotalExpectation);
		}



		cout << "rank : " << world_rank << endl;
		cout << "size : " << world_size << endl;

		TimeEnd = MPI_Wtime();
		TotalTime = TimeEnd = TimeStart;
		cout << "Time = " << TotalTime << endl;
	}
	ofile.close();  // close output file

	MPI_Finalize();

	ofstream FS;
	string fname;
	fname = "P(E)";
	fname += to_string(InitialTemp);
	fname += "_";
	fname += fileout;
	FS.open(fname);
	FS << " P(E) : " << endl;
	int size = E_vector.size();
	double ex1 = 0;
	double ex2 = 0;
	double pr, var;
	for (int c = 0; c < size; c++) {
		pr = E_vector[c].count / (double)(MCcycles);
		FS << E_vector[c].eng << ",       " << pr << endl;
		ex1 += E_vector[c].eng * pr;
		ex2 += E_vector[c].eng * E_vector[c].eng * pr;
	}

	FS << "variance : " << (ex2 - ex1*ex1) / (NSpins*NSpins) << endl;
	FS.close();
	
	/*
	int size = E_vector.size();
	double ex1 = 0;
	double ex2 = 0;
	double pr, var;
	for (int c = 0; c < size; c++) {
		pr = E_vector[c].count / (double)(MCcycles);
		cout << E_vector[c].eng << ",       " << pr << endl;
		ex1 += E_vector[c].eng * pr;
		ex2 += E_vector[c].eng * E_vector[c].eng * pr;
	}

	cout << "variance : " << ex2 - ex1*ex1 << endl;
	*/
	return 0;
}



// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, int MCcycles, double Temperature, double* ExpectationValues, long& idum, int myloop_begin, int myloop_end)
{
	/*
	// Initialize the seed and call the Mersienne algo
	std::random_device rd;
	std::mt19937_64 gen(rd());
	// Set up the uniform distribution for x \in [[0, 1]
	std::uniform_real_distribution<double> RandomNumberGenerator(0.0, 1.0);
	*/

	// Initialize the lattice spin values
	double** SpinMatrix = new double*[NSpins];
	for (int i = 0; i < NSpins; i++) {
		SpinMatrix[i] = new double[NSpins];
		for (int j = 0; j < NSpins; j++) {
			SpinMatrix[i][j] = 0;
		}
	}
	//    initialize energy and magnetization 
	double Energy = 0.;     double MagneticMoment = 0.;
	// initialize array for expectation values
	InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment);
	// setup array for possible energy changes
	double* EnergyDifference = new double[17];
	for (int i = 0; i < 17; i++) {
		EnergyDifference[i] = 0;
	}
	for (int de = -8; de <= 8; de += 4) EnergyDifference[de + 8] = exp(-de / Temperature);
	double config[2] = { -1.0, 1.0 };
	int count_config = 0;
	// Start Monte Carlo cycles
	for (int cycles = myloop_begin; cycles <= myloop_end; cycles++) {
		// The sweep over the lattice, looping over all spin sites
		for (int x = 0; x < NSpins; x++) {
			for (int y = 0; y < NSpins; y++) {
				int ix = (int)(ran2(&idum)*(double)NSpins);
				int iy = (int)(ran2(&idum)*(double)NSpins);
				int deltaE = 2 * SpinMatrix[ix][iy] * (SpinMatrix[ix][PeriodicBoundary(iy, NSpins, -1)] + SpinMatrix[PeriodicBoundary(ix, NSpins, -1)][iy] + SpinMatrix[ix][PeriodicBoundary(iy, NSpins, 1)] + SpinMatrix[PeriodicBoundary(ix, NSpins, 1)][iy]);
				if (ran2(&idum) <= EnergyDifference[deltaE + 8]) {
					SpinMatrix[ix][iy] *= -1.0;  // flip one spin and accept new spin config
					MagneticMoment += (double)2 * SpinMatrix[ix][iy];
					Energy += (double)deltaE;
					count_config+=1;
				}
			}

			//count to get probability
			if (cycles >= MCcycles) {
				if (cycles == MCcycles) {
					e_count e;
					e.eng = Energy;
					e.count++;
					E_vector.push_back(e);
				}
				else {
					int size = E_vector.size();
					for (int c = 0; c < size; c++) {
						if (Energy == E_vector[c].eng) {
							E_vector[c].count++;
						}
						else {
							e_count ee;
							ee.eng = Energy;
							ee.count++;
							E_vector.push_back(ee);
						}
					}
				}
			}
		}
		

		// update expectation values  for local node
		ExpectationValues[0] += Energy;
		ExpectationValues[1] += Energy*Energy;
		ExpectationValues[2] += MagneticMoment;
		ExpectationValues[3] += MagneticMoment*MagneticMoment;
		ExpectationValues[4] += fabs(MagneticMoment);
		ExpectationValues[5] += count_config;


	}
	cout << "count_config" <<  count_config << endl;
} // end of Metropolis sampling over spins

  // function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, double** SpinMatrix, double& Energy, double& MagneticMoment)
{
	/*
	// setup spin matrix and initial magnetization
	for (int x = 0; x < NSpins; x++) {
		for (int y = 0; y < NSpins; y++) {
			SpinMatrix[x][y] = 1.0; // spin orientation for the ground state
			MagneticMoment += (double)SpinMatrix[x][y];
		}
	}
	*/
	
	random_device seed;
	mt19937 engine(seed());
	uniform_real_distribution<double> values(0.0, 1.0);
	for (int i = 0; i < NSpins; ++i) {
		for (int j = 0; j < NSpins; ++j)
		{
			if (values(engine) < 0.5) { SpinMatrix[i][j] = 1; }
			else { SpinMatrix[i][j] = -1; }
		}
	}
	
	// setup initial energy
	for (int x = 0; x < NSpins; x++) {
		for (int y = 0; y < NSpins; y++) {
			Energy -= (double)SpinMatrix[x][y] *
				(SpinMatrix[PeriodicBoundary(x, NSpins, -1)][y] +
					SpinMatrix[x][PeriodicBoundary(y, NSpins, -1)]);
		}
	}
}// end function initialise



void WriteResultstoFile(int NSpins, int MCcycles, double temperature, double* ExpectationValues)
{
	double norm = 1.0 / ((double)(MCcycles));  // divided by  number of cycles 
	double E_ExpectationValues = ExpectationValues[0] * norm;
	double E2_ExpectationValues = ExpectationValues[1] * norm;
	double M_ExpectationValues = ExpectationValues[2] * norm;
	double M2_ExpectationValues = ExpectationValues[3] * norm;
	double Mabs_ExpectationValues = ExpectationValues[4] * norm;
	// all expectation values are per spin, divide by 1/NSpins/NSpins
	double Evariance = (E2_ExpectationValues - E_ExpectationValues*E_ExpectationValues) / NSpins / NSpins;
	double Mvariance = (M2_ExpectationValues - Mabs_ExpectationValues*Mabs_ExpectationValues) / NSpins / NSpins;
	ofile << setiosflags(ios::showpoint | ios::uppercase);
	ofile << setw(15) << setprecision(8) << temperature;
	ofile << setw(15) << setprecision(8) << E_ExpectationValues / NSpins / NSpins;
	ofile << setw(15) << setprecision(8) << Evariance / temperature / temperature;
	ofile << setw(15) << setprecision(8) << M_ExpectationValues / NSpins / NSpins;
	ofile << setw(15) << setprecision(8) << Mvariance / temperature;
	ofile << setw(15) << setprecision(8) << Mabs_ExpectationValues / NSpins / NSpins << endl;
	ofile << setw(15) << setprecision(8) << ExpectationValues[5] << endl;//configuration counts

	cout << "E_variance : " << Evariance << endl;

	/*
	//compare probability and cumputed E variance
	double e1 = 0;
	double e2 = 0;
	for (int c = 0; c <  E_vector.size(); c++){
		e1 += E_vector[c].eng * E_vector[c].count / 0.5*MCcycles;
		e2 += E_vector[c].eng * E_vector[c].eng * E_vector[c].count / 0.5*MCcycles;
		cout << e1 << endl;
		system("pause");
	}
	cout << "computed Energy variance : " << Evariance << endl;
	cout << "calculated Energy variance : " << e2 - e1*e1 << endl;
	*/
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
	int            j;
	long           k;
	static long    idum2 = 123456789;
	static long    iy = 0;
	static long    iv[NTAB];
	double         temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum = 1;
		else             *idum = -(*idum);
		idum2 = (*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ1;
			*idum = IA1*(*idum - k*IQ1) - k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB)  iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ1;
	*idum = IA1*(*idum - k*IQ1) - k*IR1;
	if (*idum < 0) *idum += IM1;
	k = idum2 / IQ2;
	idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j = iy / NDIV;
	iy = iv[j] - idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp = AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX