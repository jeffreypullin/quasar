/*
 * MIT License
 *
 * Copyright (c) 2024 Jeffrey Pullin
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
        return;
    }

    std::string line;
    std::vector<std::string> tokens;
    size_t index = 0;
    while (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, "\t ");
        if (tokens.size() != 6) {
            std::cerr << "Error: Invalid BIM file format." << std::endl;
            return;
        }
        chrom.push_back(std::stoi(tokens[0]));
        id.push_back(tokens[1]);
        pos.push_back(std::stoul(tokens[3]));
        allele1.push_back(tokens[4]);
        allele2.push_back(tokens[5]);
        this->index.push_back(index);
        index++;
    }

    n_snps = index;
    std::cout << "Number of SNPs read: " << n_snps << std::endl;
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
        return;
    }

    std::string line;
    std::vector<std::string> tokens;
    while (std::getline(file, line)) {
        tokens = string_split(line, "\t ");
        if (tokens.size() < 6) {
            std::cerr << "Error: Invalid FAM file format." << std::endl;
            return;
        }

        sample_ids.push_back(tokens[1]);
    }
    n_samples = sample_ids.size();

    std::cout << "Number of samples read: " << n_samples << std::endl;

    file.close();
}

void GenoData::prepare_bed_file() {
    std::string bed_file = bed_prefix + ".bed";
    std::ifstream file(bed_file, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open BED file " << bed_file << std::endl;
        return;
    }

    char magic[3];
    file.read(magic, 3);
    if (magic[0] != 0x6C || magic[1] != 0x1B) {
        std::cerr << "Error: Invalid BED file magic number." << std::endl;
        return;
    }
    if (magic[2] != 0x01) {
        std::cerr << "Error: BED file is not in SNP-major mode." << std::endl;
        return;
    }

    file.close();
}

void GenoData::read_bed_file() {
    std::string bed_file = bed_prefix + ".bed";
    std::ifstream file(bed_file, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open BED file " << bed_file << std::endl;
        return;
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
                    case 0: genotype_matrix(i, j + k) = 2; break;
                    case 1: genotype_matrix(i, j + k) = -1; break;
                    case 2: genotype_matrix(i, j + k) = 1; break;
                    case 3: genotype_matrix(i, j + k) = 0; break;
                } 
            }
        }
    }

    // Transpose to samples x snps.
    genotype_matrix.transposeInPlace();

    std::cout << "Genotype matrix read successfully." << std::endl;

    file.close();
}

void GenoData::slice_samples(std::vector<std::string>& sample_ids) {
    
    Eigen::VectorXi rows;
    rows.resize(sample_ids.size());
    for (int i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            rows(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in phenotype data." << std::endl;
            return;
        }
    }
    // Create a new matrix with the selected columns
    Eigen::MatrixXd sliced_genotype_matrix(rows.size(), genotype_matrix.cols());
    for (size_t i = 0; i < rows.size(); ++i) {
        sliced_genotype_matrix.row(i) = genotype_matrix.row(rows[i]);
    }
    
    this->genotype_matrix = sliced_genotype_matrix;
    this->sample_ids = sample_ids;
    n_samples = sample_ids.size();
}   

void GenoData::standardise(Eigen::MatrixXd& X) {

    std::cout << "Projecting out covariates..." << std::endl;
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(X);
    Eigen::MatrixXd Q = qr.householderQ() * Eigen::MatrixXd::Identity(X.rows(), X.cols());

    Eigen::MatrixXd QTG = Q.transpose() * genotype_matrix;
    genotype_matrix.noalias() -= Q * QTG;

    std::cout << "Standardising..." << std::endl;
    for (int i = 0; i < genotype_matrix.cols(); ++i) {
        double mean = genotype_matrix.col(i).mean();
        genotype_matrix.col(i).array() -= mean;

        double stddev = std::sqrt(genotype_matrix.col(i).array().square().sum() / (n_samples - 1));
        if (stddev == 0) { 
            stddev = 1;
        }
        genotype_matrix.col(i) /= stddev;
    }
}