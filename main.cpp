/*
Calculating mean energy, mean magnetization,
specific heat and susceptibility as function of temperature with PBC
through Ising model in 2-dimensions, using Metropolis algorithm
*/
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace std;

ofstream ofile;

//inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add)
{
    return (i + limit + add) % (limit);
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



int main (int argc, char* argv[]){ //main program

if (argc <= 1){
    cout << "Naah. Not enough input: " << argv[0] << ". Include output name on same line." << endl;
    exit(1);
    } //end if
else{
    char* outfilename = argv[1];
    }  //end else

//Set up array with possible energy-changes
double w[17], temperature;
for (int de = -8; de <= 8; de++) w[de+8] = 0;   //index 0 corresponds to deltaE = -8
for (int de = -8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);   //deltaE can only in-/decrease by 4

//Get matrix with spins

//Set initial energy and magnetization as well as random starting point

//Calculate energy and magnetization

//Flip one spin



} //end of main
