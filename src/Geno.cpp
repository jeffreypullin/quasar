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

void read_bim_file(Param* params, std::vector<SNP>* snps_info) {

    std::string bim_file = params->bed_prefix + ".bim";
    std::ifstream file(bim_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open BIM file " << bim_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;
    while (std::getline(file, line)) {
        tokens = string_split(line, "\t");
        if (tokens.size() != 6) {
            std::cerr << "Error: Invalid BIM file format." << std::endl;
            return;
        }

        SNP snp;
        snp.chrom = std::stoi(tokens[0]);
        snp.id = tokens[1];
        snp.pos = std::stoul(tokens[3]);
        snp.allele1 = tokens[4];
        snp.allele2 = tokens[5];
        // Placeholder, MAF calculation not included in this function.
        snp.MAF = 0.0;

        snps_info->push_back(snp);
    }

    std::cout << "Number of SNPs read: " << snps_info->size() << std::endl;

    file.close();

}

void read_fam_file(Param* params, std::vector<std::string>* sample_ids, int* n_samples) {

    std::string fam_file = params->bed_prefix + ".fam";
    std::ifstream file(fam_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open FAM file " << fam_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;
    while (std::getline(file, line)) {
        tokens = string_split(line, "\t");
        if (tokens.size() < 6) {
            std::cerr << "Error: Invalid FAM file format." << std::endl;
            return;
        }

        sample_ids->push_back(tokens[0] + "_" + tokens[1]);
    }
    *n_samples = sample_ids->size();

    std::cout << "Number of samples read: " << sample_ids->size() << std::endl;

    file.close();
}

void prepare_bed_file(Param* params) {
    std::string bed_file = params->bed_prefix + ".bed";
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

void read_bed_file_chunk(Param* params, Eigen::MatrixXd* genotype_matrix, int* n_samples, std::vector<size_t>& snp_indices) {
    std::string bed_file = params->bed_prefix + ".bed";
    std::ifstream file(bed_file, std::ios::binary);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open BED file " << bed_file << std::endl;
        return;
    }

    // Skip the magic number and mode.
    file.seekg(3, std::ios::beg);

    int n = *n_samples;
    genotype_matrix->resize(snp_indices.size(), n);

    for (size_t i = 0; i < snp_indices.size(); ++i) {
        size_t snp_index = snp_indices[i];
            file.seekg(3 + snp_index * ((n + 3) / 4), std::ios::beg);

        for (size_t j = 0; j < n; j += 4) {
            char byte;
            file.read(&byte, 1);

            for (size_t k = 0; k < 4 && j + k < n; ++k) {
                char genotype = (byte >> (2 * k)) & 0x03;
                switch (genotype) {
                    // Homozygous for allele1.
                    case 0x00:
                        (*genotype_matrix)(i, j + k) = 0;
                        break;
                    // Heterozygous.
                    case 0x01:
                        (*genotype_matrix)(i, j + k) = 1;
                        break;
                    // Homozygous for allele2.
                    case 0x02:
                        (*genotype_matrix)(i, j + k) = 2;
                        break;
                    // Missing genotype.
                    case 0x03:
                        (*genotype_matrix)(i, j + k) = -1;
                        break;
                }
            }
        }
    }

    std::cout << "Genotype matrix read successfully." << std::endl;

    file.close();
}
