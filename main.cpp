/*
Calculating mean energy, mean magnetization,
specific heat and susceptibility as function of temperature with PBC
through Ising model in 2-dimensions, using Metropolis algorithm

#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
ofstream ofile;

//function to get values from command line
void ReadInput(double initialTemperature, double finalTemperature, double temp_step, int cycles, int n_spins){
    cout << "Initial temperature: " << endl;
    cin >> initialTemperature;
    cout << "Final temperature: " << endl;
    cin >> finalTemperature;
    cout << "Temperature increments by: " << endl;
    cin >> temp_step;
    cout << "Number of MC-cycles: " << endl;
    cin >> cycles;
    cout << "Number of spins in lattice: " << endl;
    cin >> n_spins;
}

//make initial spin-lattice (grid) and initial magnetization
void Initialize(int n_spins, double temperature, int **spin_matrix, double& E, double& M)
{
    for (int y = 0; y < n_spins; y++) {
        for (int x = 0; x < n_spins; x++) {
            if (temperature < 1.5) {
                spin_matrix[y][x] = 1;  //spin orientation for ground state
                M += (double) spin_matrix[y][x];
            } //end if-test
        }
    }   //end for-loops
// How is the spin-matrix for T>1.5 set up??
    //set up initial energy
    for (int y = 0; y < n_spins; y++) {
        for (int x = 0; x < n_spins; x++) {
            E -= (double) spin_matrix[y][x]*
                    (spin_matrix[periodic(y, n_spins, -1)][x] +
                     spin_matrix[y][periodic(x, n_spins, -1)]);
        }
    }
} //end initialize

//Metropolis algorithm
void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double& M, double *w){
    for (int x = 0; x <= n_spins; x++){
        for (int y = 0; y <= n_spins; y++){
        //finding a random position for the spin-flip
        int ix = (int) (ran1(&idum)*(double)n_spins);
        int iy = (int) (ran1(&idum)*(double)n_spins);
        int deltaE = 2*spin_matrix[iy][ix]*(spin_matrix[iy][ix-1] +
                                        spin_matrix[iy][ix+1] +
                                        spin_matrix[iy+1][ix] +
                                        spin_matrix[iy-1][ix]);
        if (deltaE > 0){
            if (ran1(&idum) < deltaE){
                spin_matrix[iy][ix] *= -1;  //flip the spin
                M += (double) 2*spin_matrix[iy][ix];
                E += deltaE;
                }
            else{
                spin_matrix[iy][ix] += 0;   //no change
                }
            }
        else{   //if deltaE <= 0
            spin_matrix[iy][ix] += 0;
            }
        }
    }   //end for-loops
} //end Metropolis


//main program
int main (int argc, char* argv[]){

int n_spins = 2;
int cycles = 10;
double initialTemperature = 1.0;
double finalTemperature = 2.0;
double temp_step = 0.1;
int spin_matrix[n_spins][n_spins];
for (int i = 0; i < n_spins; i++){
    for (int j = 0; j<n_spins; j++){
    spin_matrix[j][i] = 1;
    cout << spin_matrix[j][i] << endl;
    }
}
if (argc <= 1){
    cout << "Naah. Not enough input: " << argv[0] << ". Include output name on same line." << endl;
    exit(1);
    } //end if
else{
    char* outfilename = argv[1];
    }  //end else

//get initial values from command line (# cycles, temperature etc.)
//ofile.open(outfilename);

ReadInput(initialTemperature, finalTemperature, temp_step, cycles, n_spins);

for (double temperature = 1.0; temperature <= 2.0; temperature += 0.1){
    int E = 0;  //initial energy and magnetization
    int M = 0;

    //Set up array with possible energy-changes
    double w[17];
    for (int de = -8; de <= 8; de++) w[de+8] = 0;   //index 0 corresponds to deltaE = -8
    for (int de = -8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);   //deltaE can only in-/decrease by 4

    //initialise array for expectation values
    double average[5];
    for (int i = 0; i < 5; i++) average[i] = 0;
//    Initialize(n_spins, temperature, spin_matrix, E, M);
}
//Get matrix with spins

//Set initial energy and magnetization as well as random starting point

//Calculate energy and magnetization

//Flip one spin


return 0;
} //end of main
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

ofstream ofile;


int** setup_matrix(int L) {
    /*This function will set up an initial (L+2)x(L+2)-lattice where
    every element equals 1. */
    int** mat = new int*[L+2];

    for (int i = 0; i < L+2; i++){
        mat[i] =new int[L];
    }
    for (int i = 0; i < L+2; i++){
        for (int j = 0; j < L+2; j++){
            mat[j][i] = 1; //set every element = 1

        }
    } //end for-loops
    return mat;
}

void metropolis(long &E, long &M, int** spin_matrix, int L, double beta) {

    //find random element
    for (int x = 1; x < L+1; x++){
        for (int y = 1; y < L+1; y++){
            int initx = (int) rand() % L+1;
            int inity = (int) rand() % L+1;

            //calculate change in energy
            int deltaE = (2*spin_matrix[inity][initx])*(spin_matrix[inity-1][initx] + spin_matrix[inity+1][initx] +
                                                        spin_matrix[inity][initx-1] + spin_matrix[inity][initx+1]);
            bool accepted = false;
            /*if energy is higher than last state, decide if new state is accepted
            through probability. The probability for acceptance decreases as the change
            in energy increases.*/
            if (deltaE > 0){
                double w = exp(-beta/deltaE);

                double r = rand() / (double)RAND_MAX;
                if (r < w){
                    accepted = true;
                }
            } else {
                accepted = true;
            }

            if(accepted) {
                spin_matrix[inity][initx] *= -1;
                E += deltaE;
                M += 2*(spin_matrix[inity][initx]);

                if (int initx = 1){spin_matrix[L][inity] = spin_matrix[initx][inity]; }
                if (int initx = L){spin_matrix[0][inity] = spin_matrix[initx][inity]; }
                if (int inity = 1){spin_matrix[initx][L] = spin_matrix[initx][inity]; }
                if (int inity = L){spin_matrix[initx][0] = spin_matrix[initx][inity]; }
            }
        }
    }
}

int calculateEnergy(int **spin_matrix, int L) {
    long E = 0;
    for (int i = 1; i < L+1; i++){
        for (int j = 1; j < L+1; j++){
            E -= spin_matrix[j][i]*(spin_matrix[j-1][i] + spin_matrix[j][i-1]);
        }
    }
    return E;
}

int calculateMagnetization(int **spin_matrix, int L) {
    int M = 0;
    for (int i = 1; i < L+1; i++){
        for (int j = 1; j < L+1; j++){
            M += spin_matrix[j][i];
        }
    }
    return M;
}

void runMC(int **spin_matrix, int L, int NCycles, double &ESum, double &ESquaredSum, double &MSum, double &MSquaredSum, double beta) {

    long E = calculateEnergy(spin_matrix, L);
    long M = calculateMagnetization(spin_matrix, L);
    ESum = 0;
    ESquaredSum = 0;
    MSum = 0;
    MSquaredSum = 0;

    for (double cycle = 1; cycle <= NCycles; cycle++){
        metropolis(E, M, spin_matrix, L, beta);
        ESum += E;
        ESquaredSum += E*E;
        MSum += M;
        MSquaredSum += M*M;
    }
    ESum /= NCycles;
    ESquaredSum /= NCycles;
    MSum /= NCycles;
    MSquaredSum /= NCycles;

}

void output(double temperature, double NCycles, double cv, double chi, double ESum, double MSum) {
    //write interesting values to file
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(6) << setprecision(2) << temperature;
    ofile << setw(15) << setprecision(8) << NCycles;
    ofile << setw(15) << setprecision(8) << cv;
    ofile << setw(15) << setprecision(8) << chi;
    ofile << setw(15) << setprecision(8) << ESum;
    ofile << setw(15) << setprecision(8) << MSum;
    ofile << "\n";
}

int main() {
    double minTemperature = 1.0;
    double finalTemperature = 3.0;
    double T_step = 0.1;

    ofile.open("testings");
    ofile << setw(6) << setprecision(8) << "temp";
    ofile << setw(15) << setprecision(8) << "NCycles";
    ofile << setw(10) << setprecision(8) << "cv";
    ofile << setw(15) << setprecision(8) << "chi";
    ofile << setw(15) << setprecision(8) << "ESum";
    ofile << setw(15) << setprecision(8) << "MSum";
    ofile << "\n";

    for (double temperature = minTemperature; temperature <= finalTemperature; temperature += T_step){
        double beta = 1/temperature;
        int L = 2;
        int seed = -1;

        //initialize random seed
        srand(seed);

        //call function setup_matrix
        int** spin_matrix = setup_matrix(L);

        double NCycles = 1e5;
        double meanEnergy = 0;
        double meanEnergySquared = 0;
        double meanM = 0;
        double meanMSquared = 0;

        // Thermalize (reach equilibrium)
        runMC(spin_matrix, L, 1e4, meanEnergy, meanEnergySquared, meanM, meanMSquared, beta);

        // Start sampling for reals
        runMC(spin_matrix, L, 1e5, meanEnergy, meanEnergySquared, meanM, meanMSquared, beta);

        double cv = (beta/temperature) * (meanEnergySquared - meanEnergy*meanEnergy) / (L*L);
        double chi = beta * (meanMSquared - meanM*meanM);

        output(temperature, NCycles, cv, chi, meanEnergy, meanM);

        cout << temperature << endl;
    }
    return 0;
}
