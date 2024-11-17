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
        tokens = string_split(line, " ");
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
        tokens = string_split(line, " ");
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
        size_t snp_index = index[i];

        file.seekg(3 + snp_index * ((n_samples + 3) / 4), std::ios::beg);

        for (size_t j = 0; j < n_samples; j += 4) {
            char byte;
            file.read(&byte, 1);

            for (size_t k = 0; k < 4 && j + k < n_samples; ++k) {
                char genotype = (byte >> (2 * k)) & 0x03;
                switch (genotype) {
                    // Homozygous for allele1.
                    case 0x00:
                        genotype_matrix(i, j + k) = 0;
                        break;
                    // Heterozygous.
                    case 0x01:
                        genotype_matrix(i, j + k) = 1;
                        break;
                    // Homozygous for allele2.
                    case 0x02:
                        genotype_matrix(i, j + k) = 2;
                        break;
                    // Missing genotype.
                    case 0x03:
                        genotype_matrix(i, j + k) = -1;
                        break;
                }
            }
        }
    }

    std::cout << "Genotype matrix read successfully." << std::endl;

    file.close();
}

void GenoData::slice_samples(std::vector<std::string>& sample_ids) {
    
    // Create a vector to store the indices of the samples we want to keep
    Eigen::VectorXi cols;
    cols.resize(sample_ids.size());
    for (int i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            cols(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in phenotype data." << std::endl;
            return;
        }
    }
    // Create a new matrix with the selected columns
    Eigen::MatrixXd sliced_genotype_matrix(genotype_matrix.rows(), cols.size());
    for (size_t i = 0; i < cols.size(); ++i) {
        sliced_genotype_matrix.col(i) = genotype_matrix.col(cols[i]);
    }
    
    this->genotype_matrix = sliced_genotype_matrix;
    this->sample_ids = sample_ids;
    n_samples = sample_ids.size();
}   

void GenoData::compute_variant_var() {
    var.resize(n_snps);
    for (size_t i = 0; i < n_snps; ++i) {
        double mean = genotype_matrix.row(i).mean();
        var[i] = (genotype_matrix.row(i).array().pow(2).sum() - (mean * mean) * n_samples) / (n_samples - 1);
    }
}