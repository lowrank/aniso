/*
 * output csv subroutine.
 *
 * copyright@ Yimin Zhong. yzhong@math.utexas.edu.
 */

#ifndef MATLAB_IO_H
#define MATLAB_IO_H

#include "bbfmm/utils.h"
#include "bbfmm/bbfmm.h"
#include "bbfmm/linalg.h"

void write_to_csv(vector<scalar_t> &data, std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    for (auto it : data) {
        outfile << it << "\n";
    }
    outfile.close();
}


void write_to_csv(bbfmm::Vector &data, std::string filename) {
    std::ofstream outfile;
    outfile.open(filename);
    for (int id = 0; id < data.row(); ++id) {
        outfile << std::scientific << std::setprecision(32) << std::setw(32) << data(id) << "\n";
    }
    outfile.close();
}

void write_to_csv(vector<point> &data, std::string filename, std::string sep = " ") {
    std::ofstream outfile;
    outfile.open(filename);
    for (auto it : data) {
        outfile << std::scientific << std::setprecision(32) << std::setw(32) << it.x << sep << std::setw(32) << it.y
                << "\n";
    }
    outfile.close();
}

void load_csv(bbfmm::Vector &data, std::string filename) {
    std::ifstream infile;
    infile.open(filename);
    if (!infile.is_open()) {
        std::cout << "File not open." << std::endl;
    }

    int row = 0;
    double num = 0.;
    while (infile >> num) {
        data(row) = num;
        row++;
    }
}


#endif //MATLAB_IO_H
