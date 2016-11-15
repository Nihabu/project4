#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <armadillo>
#include <ctime>
using namespace arma;
using namespace std;

ofstream ofile;

int** setup_matrix(int L) {
    /*This function will set up an initial (L+2)x(L+2)-lattice where
    every element equals 1. */
    int** mat = new int*[L+2];

    for (int i = 0; i < L+2; i++){
        mat[i] =new int[L+2];
    }
    for (int i = 0; i < L+2; i++){
        for (int j = 0; j < L+2; j++){
            //mat[j][i] = 1; //set every element = 1
            //following lines sets up a random matrix

            double rand_int = rand()/(double)RAND_MAX;
            if (rand_int <= 0.5){ mat[j][i] = -1; }
            else { mat[j][i] = 1; }

        }
    } //end for-loops
    return mat;
}

void metropolis(double &E, double &M, int** spin_matrix, int L, vec exponentials, double &totAccepted) {
    totAccepted = 0;
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
                double w = exponentials(deltaE+8);
                double r = rand() / (double)RAND_MAX;
                if (r < w){
                    accepted = true;
                }
            } else {                
                accepted = true;

            }
            if(accepted) {
                totAccepted += 1;
                spin_matrix[inity][initx] *= -1;
                //cout << "x: " << initx << "   y: " << inity << endl;
                E += deltaE;
                M += 2*spin_matrix[inity][initx];

                if (int initx = 1){spin_matrix[inity][L+1] = spin_matrix[inity][initx]; }
                if (int initx = L){spin_matrix[inity][0] = spin_matrix[inity][initx]; }
                if (int inity = 1){spin_matrix[L+1][initx] = spin_matrix[inity][initx]; }
                if (int inity = L){spin_matrix[0][initx] = spin_matrix[inity][initx]; }
            }
        }
    }
}

int calculateEnergy(int **spin_matrix, int L) {
    double E = 0;
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

void output(double temperature, double NCycles, double cv, double chi, double meanEnergy, double meanMAbs, double totAccepted, double meanESquared) {
    //write interesting values to file
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(6) << setprecision(3) << NCycles;
    ofile << setw(15) << setprecision(8) << temperature;
    ofile << setw(15) << setprecision(8) << cv;
    ofile << setw(15) << setprecision(8) << chi;
    ofile << setw(15) << setprecision(8) << meanEnergy;
    ofile << setw(15) << setprecision(8) << meanMAbs;
    ofile << setw(15) << setprecision(8) << totAccepted;
    ofile << setw(15) << setprecision(8) << meanESquared;
    ofile << "\n";
}

void runMC(int **spin_matrix, int L, int NCycles, double &ESum, double &ESquaredSum, double &MSum, double &MSquaredSum, double &MAbs, vec &exponentials, double temperature, double totAccepted) {
    double E = calculateEnergy(spin_matrix, L);
    double M = calculateMagnetization(spin_matrix, L);
    ESum = 0;
    ESquaredSum = 0;
    MSum = 0;
    MSquaredSum = 0;
    MAbs = 0;

    for (double cycle = 1; cycle <= NCycles; cycle++){
        metropolis(E, M, spin_matrix, L, exponentials, totAccepted);
        ESum += E;
        ESquaredSum += E*E;
        MSum += M;
        MSquaredSum += M*M;
        MAbs += fabs(M);

    }
    ESum /= NCycles;
    ESquaredSum /= NCycles;
    MSum /= NCycles;
    MSquaredSum /= NCycles;
    MAbs /= NCycles;
}

int main() {
    double minTemperature = 2.2;
    double finalTemperature = 2.3;
    double T_step = 0.01;
    //double temperature = 2.4;

    ofile.open("filename");

    for (double temperature = minTemperature; temperature <= finalTemperature+T_step; temperature += T_step){
        cout << temperature << endl;
        double beta = 1.0/temperature;
        int L = 40;
        //int seed = -2;
        //initialize random seed
        srand(time(NULL));

        //make code more efficient by pre-calculating exponentials
        vec exponentials;
        exponentials.zeros(17);
        exponentials(-8+8) = exp(-beta*(-8));
        exponentials(-4+8) = exp(-beta*(-4));
        exponentials(0+8) = 1;
        exponentials(4+8) = exp(-beta*4);
        exponentials(8+8) = exp(-beta*8);

        //call function setup_matrix
        int** spin_matrix = setup_matrix(L);

        double meanEnergy = 0;
        double meanEnergySquared = 0;
        double meanM = 0;
        double meanMSquared = 0;
        double meanMAbs = 0;
        double totAccepted = 0;
        double NCycles = 1e6;

        // Thermalize (reach equilibrium)
        runMC(spin_matrix, L, 10000, meanEnergy, meanEnergySquared, meanM, meanMSquared, meanMAbs, exponentials, temperature, totAccepted);

        // Start sampling for reals
        runMC(spin_matrix, L, NCycles, meanEnergy, meanEnergySquared, meanM, meanMSquared, meanMAbs, exponentials, temperature, totAccepted);

        double cv = (beta/temperature) * (meanEnergySquared - meanEnergy*meanEnergy) / (L*L);
        double chi = beta  * (meanMSquared - meanMAbs*meanMAbs) / (L*L);
        double EperSpin = meanEnergy/(L*L);
        double MperSpin = meanMAbs/(L*L);

        output(temperature, NCycles, cv, chi, EperSpin, MperSpin, totAccepted, meanEnergySquared);
    }
    return 0;
}
