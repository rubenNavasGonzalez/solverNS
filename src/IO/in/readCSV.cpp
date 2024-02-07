//
// Created by ruben on 6/02/24.
//

#include <fstream>
#include <sstream>
#include "readCSV.h"


std::vector< std::vector<double> > readCSV(std::string filename) {


    // Auxiliary variables declaration
    std::string line;


    // .csv data matrix initialization
    std::vector< std::vector<double> > data;


    // Create the inFile
    std::ifstream inFile(filename);


    // Check if file is opened
    if (inFile.is_open()) {


        // Read the header (skip it)
        std::getline(inFile, line);


        // Read the numerical data row by row
        while ( std::getline(inFile, line) ) {

            std::vector<double> row;
            std::stringstream linestream(line);
            std::string value;


            // Read column by column
            while ( std::getline(linestream, value, ',') ) {

                row.push_back(std::stod(value));
            }

            data.push_back(row);
        }

    } else {

        printf("Error. Unable to open the .csv file \n");
    }


    return data;
}