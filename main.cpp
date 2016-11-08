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
//#include "lib.h"
using namespace std;

//metropolis algorithm

int main() {
    //set up initial spin_matrix
    int L = 6;  //dimensions
    double temperature = 1.5; //temperature
    int idum = -1;

    int** spin_matrix = new int*[L];
    for (int x = 0; x < L; x++){
        spin_matrix[x] = new int[L];
    }

    for (int i = 0; i < L; i++){
        for (int j = 0; j < L; j++){
            spin_matrix[j][i] = 1;
        }
    }//end for-loops

    /* initialize random seed: */
    srand (idum);
    //find random element
    for (int x = 0; x < L; x++){
        for (int y = 0; y < L; y++){
            int initx = (int) rand() % L;
            int inity = (int) rand() % L;
            int deltaE = 2*spin_matrix[inity][initx]*(spin_matrix[inity-1][initx] + spin_matrix[inity+1][initx] +
                                                        spin_matrix[inity][initx-1] + spin_matrix[inity][initx+1]);
            if (deltaE > 0){
                if (rand() < deltaE){
                    spin_matrix[inity][initx] *= -1;
                }
                else{
                    //calculate what happens if the test is not a success
                }
            }
            else{
                spin_matrix[inity][initx] *= -1;

            }
        }
    }

    //array with possible energy changes (8, 4, 0, -4, -8)
    double w[17]; //array with size 17 (-8 to 8)
    for (int de = -8; de <= 8; de++) w[de+8] = 0;   //index 0 corresponds to energychange -8
    for (int de = -8; de <= 8; de += 4) w[de+8] = exp(-de/temperature); //deltaE can only in-/decrease by 4
    cout << w[1] << endl;

    return 0;
}
