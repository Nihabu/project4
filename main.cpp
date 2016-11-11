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

void output(double temperature, double NCycles, double cv, double chi, double meanEnergy, double meanM) {
    //write interesting values to file
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(6) << setprecision(3) << temperature;
    ofile << setw(15) << setprecision(8) << NCycles;
    ofile << setw(15) << setprecision(8) << cv;
    ofile << setw(15) << setprecision(8) << chi;
    ofile << setw(15) << setprecision(8) << meanEnergy;
    ofile << setw(15) << setprecision(8) << meanM;
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

    for (double temperature = minTemperature; temperature <= finalTemperature+T_step; temperature += T_step){
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
        runMC(spin_matrix, L, 1e6, meanEnergy, meanEnergySquared, meanM, meanMSquared, beta);

        double cv = (beta/temperature) * (meanEnergySquared - meanEnergy*meanEnergy) / (L*L);
        double chi = beta * (meanMSquared - meanM*meanM);

        output(temperature, NCycles, cv, chi, meanEnergy, meanM);
    }
    return 0;
}
