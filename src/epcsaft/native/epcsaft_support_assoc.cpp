#include "epcsaft_core_internal.h"

using namespace thermo_detail;

vector<double> association_site_fractions_cpp(vector<double> XA_guess, vector<double> delta_ij, double den,
    vector<double> x) {
    /**Iterate over this function in order to solve for XA*/
    int num_sites = static_cast<int>(XA_guess.size());
    vector<double> XA = XA_guess;

    int idxij = -1; // index for delta_ij
    for (int i = 0; i < num_sites; i++) {
        double summ = 0.;
        for (int j = 0; j < num_sites; j++) {
            idxij += 1;
            summ += den*x[j]*XA_guess[j]*delta_ij[idxij];
        }
        XA[i] = 1./(1.+summ);
    }

    return XA;
}


vector<double> association_site_fraction_dt_cpp(vector<double> delta_ij, double den,
    vector<double> XA, vector<double> ddelta_dt, vector<double> x) {
    /**Solve for the derivative of XA with respect to temperature.*/
    int num_sites = static_cast<int>(XA.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_sites, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites, num_sites);

    double summ;
    int ij = 0;
    for (int i = 0; i < num_sites; i++) {
        summ = 0;
        for (int j = 0; j < num_sites; j++) {
            B(i) -= x[j]*XA[j]*ddelta_dt[ij];
            A(i,j) = x[j]*delta_ij[ij];
            summ += x[j]*XA[j]*delta_ij[ij];
            ij += 1;
        }
        A(i,i) = pow(1+den*summ, 2.)/den;
    }

    Eigen::MatrixXd solution = A.lu().solve(B); //Solves linear system of equations
    vector<double> dXA_dt(num_sites);
    for (int i = 0; i < num_sites; i++) {
        dXA_dt[i] = solution(i);
    }
    return dXA_dt;
}


vector<double> association_site_fraction_dx_cpp(vector<int> assoc_num, vector<double> delta_ij,
    double den, vector<double> XA, vector<double> ddelta_dx, vector<double> x) {
    /**Solve for the derivative of XA with respect to composition, or actually
    rho_i (the molar density of component i, which equals x_i * rho).*/
    int num_sites = static_cast<int>(XA.size());
    int ncomp = static_cast<int>(assoc_num.size());
    Eigen::MatrixXd B(num_sites*ncomp, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites*ncomp, num_sites*ncomp);

    double sum1, sum2;
    int idx1 = 0;
    int ij = 0;
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < num_sites; j++) {
            sum1 = 0;
            for (int k = 0; k < num_sites; k++) {
                sum1 = sum1 + den*x[k]*(XA[k]*ddelta_dx[i*num_sites*num_sites + j*num_sites + k]);
                A(ij,i*num_sites+k) = XA[j]*XA[j]*den*x[k]*delta_ij[j*num_sites+k];
            }

            sum2 = 0;
            for (int l = 0; l < assoc_num[i]; l++) {
                sum2 = sum2 + XA[idx1+l]*delta_ij[idx1*num_sites+l*num_sites+j];
            }

            A(ij,ij) = A(ij,ij) + 1;
            B(ij) = -1*XA[j]*XA[j]*(sum1 + sum2);
            ij += 1;
        }
        idx1 += assoc_num[i];
    }

    Eigen::MatrixXd solution = A.lu().solve(B); //Solves linear system of equations
    vector<double> dXA_dx(num_sites*ncomp);
    for (int i = 0; i < num_sites*ncomp; i++) {
        dXA_dx[i] = solution(i);
    }
    return dXA_dx;
}
