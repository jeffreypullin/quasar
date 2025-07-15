/*
  Copyright (C) 2024-25 Jeffrey Pullin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 or 3 of the License
  (at your option).

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is available at
   http://www.r-project.org/Licenses/
*/

#include "Geno.hpp"
#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void GenoData::read_bim_file() {

    std::string bim_file = bed_prefix + ".bim";
    std::ifstream file(bim_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open BIM file " << bim_file << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string> tokens;
    size_t index = 0;
    while (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, "\t ");
        if (tokens.size() != 6) {
            std::cerr << "Error: Invalid BIM file format." << std::endl;
            exit(1);
        }
        chrom.push_back(std::stoi(tokens[0]));
        id.push_back(tokens[1]);
        pos.push_back(std::stoul(tokens[3]));
        alt.push_back(tokens[4]);
        ref.push_back(tokens[5]);
        this->index.push_back(index);
        index++;
    }

    n_snps = index;
    std::cout << "Number of SNPs: " << n_snps << std::endl;
    std::cout << "SNPs on chromosome(s): ";
    std::vector<int> unique_chrom = chrom;
    std::sort(unique_chrom.begin(), unique_chrom.end());
    unique_chrom.erase(std::unique(unique_chrom.begin(), unique_chrom.end()), unique_chrom.end());
    for (size_t i = 0; i < unique_chrom.size(); ++i) {
        std::cout << unique_chrom[i];
        if (i < unique_chrom.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;

    file.close();
}

void GenoData::read_fam_file() {

    std::string fam_file = bed_prefix + ".fam";
    std::ifstream file(fam_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open FAM file " << fam_file << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string> tokens;
    while (std::getline(file, line)) {
        tokens = string_split(line, "\t ");
        if (tokens.size() < 6) {
            std::cerr << "Error: Invalid FAM file format." << std::endl;
            exit(1);
        }

        sample_ids.push_back(tokens[1]);
    }
    n_samples = sample_ids.size();

    std::cout << "Number of samples: " << n_samples << std::endl;

    file.close();
}

void GenoData::prepare_bed_file() {
    std::string bed_file = bed_prefix + ".bed";
    std::ifstream file(bed_file, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open BED file " << bed_file << std::endl;
        exit(1);
    }

    char magic[3];
    file.read(magic, 3);
    if (magic[0] != 0x6C || magic[1] != 0x1B) {
        std::cerr << "Error: Invalid BED file magic number." << std::endl;
        exit(1);
    }
    if (magic[2] != 0x01) {
        std::cerr << "Error: BED file is not in SNP-major mode." << std::endl;
        exit(1);
    }

    file.close();
}

void GenoData::read_bed_file() {
    std::string bed_file = bed_prefix + ".bed";
    std::ifstream file(bed_file, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open BED file " << bed_file << std::endl;
        exit(1);
    }

    // Skip the magic number and mode.
    file.seekg(3, std::ios::beg);
    genotype_matrix.resize(n_snps, n_samples);

    for (size_t i = 0; i < n_snps; ++i) {

        for (size_t j = 0; j < n_samples; j += 4) {
            unsigned char byte;
            file.read(reinterpret_cast<char*>(&byte), 1); 
            
            for (size_t k = 0; k < 4 && j + k < n_samples; ++k) {
                int genotype = (byte >> (k << 1)) & 0x3;
                 switch(genotype) {
                    case 0: genotype_matrix(i, j + k) = 0; break;
                    case 1: genotype_matrix(i, j + k) = -1; break;
                    case 2: genotype_matrix(i, j + k) = 1; break;
                    case 3: genotype_matrix(i, j + k) = 2; break;
                } 
            }
        }
    }

    // Transpose to samples x snps.
    genotype_matrix.transposeInPlace();

    std::cout << "Genotype matrix read successfully." << std::endl;

    file.close();
}

void GenoData::run_mean_imputation() {

    for (int j = 0; j < this->genotype_matrix.cols(); ++j) {
        auto col = this->genotype_matrix.col(j);
        
        if ((col.array() == -1).any()) {
            
            double sum = 0;
            int count = 0;
            for (int i = 0; i < col.size(); ++i) {
                if (col(i) != -1) {
                    sum += col(i);
                    count++;
                }
            }
            double mean = sum / count;

            for (int i = 0; i < col.size(); ++i) {
                if (col(i) == -1) {
                    col(i) = mean;
                }
            }
        }
    }

    std::cout << "Mean imputation completed successfully." << std::endl;
}

void GenoData::slice_samples(std::vector<std::string>& sample_ids) {
    
    Eigen::VectorXi rows;
    rows.resize(sample_ids.size());
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            rows(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in phenotype data." << std::endl;
            exit(1);
        }
    }
    // Create a new matrix with the selected columns
    Eigen::MatrixXd sliced_genotype_matrix(rows.size(), genotype_matrix.cols());
    for (Eigen::Index i = 0; i < rows.size(); ++i) {
        sliced_genotype_matrix.row(i) = genotype_matrix.row(rows[i]);
    }
    
    this->genotype_matrix = sliced_genotype_matrix;
    this->sample_ids = sample_ids;
    n_samples = sample_ids.size();
}
