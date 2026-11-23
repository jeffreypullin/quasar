/*
  Copyright (C) 2024-26 Jeffrey Pullin

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
#include <cmath>
#include <algorithm>
#include <numeric>
#include <random>

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
        try {
            chrom.push_back(std::stoi(tokens[0]));
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Invalid chromosome detected: " << tokens[0] << std::endl;
            std::cerr << "quasar currently only supports autosomal chromosomes." << std::endl;
            std::cerr << "Data can be filtered with plink2: " << std::endl;
            std::cerr << "plink2 --bfile {data} --autosome --make-bed --out {autosome_data}" << std::endl;
            exit(1);
        }
        id.push_back(tokens[1]);
        pos.push_back(std::stoul(tokens[3]));
        alt.push_back(tokens[4]);
        ref.push_back(tokens[5]);
        snp_id.push_back(tokens[0] + ":" + tokens[3] + tokens[5] + "-" + tokens[4]);
        this->index.push_back(index);
        index++;
    }

    n_snps = index;
    std::cout << "Read " << format_with_commas(n_snps) << " SNPs from .bim file." << std::endl;
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

void GenoData::subset_region(int region_chr, int region_start, int region_end) {

    std::vector<size_t> keep;
    for (size_t i = 0; i < n_snps; ++i) {
        if (chrom[i] == region_chr && pos[i] >= region_start && pos[i] <= region_end) {
            keep.push_back(i);
        }
    }
    if (keep.empty()) {
        std::cerr << "Error: no variants found in region " << region_chr << ":"
                  << region_start << "-" << region_end << "." << std::endl;
        exit(1);
    }

    auto subset = [&keep](auto& v) {
        std::remove_reference_t<decltype(v)> out;
        out.reserve(keep.size());
        for (size_t i : keep) {
            out.push_back(v[i]);
        }
        v = std::move(out);
    };
    subset(chrom);
    subset(pos);
    subset(id);
    subset(alt);
    subset(ref);
    subset(snp_id);
    subset(index);

    // Genotypes are samples x SNPs, so keep only the selected columns.
    Eigen::MatrixXd sliced(genotype_matrix.rows(), keep.size());
    for (size_t k = 0; k < keep.size(); ++k) {
        sliced.col(k) = genotype_matrix.col(keep[k]);
    }
    genotype_matrix = std::move(sliced);
    n_snps = keep.size();

    std::cout << "Restricted to " << format_with_commas(n_snps) << " SNPs in region "
              << region_chr << ":" << region_start << "-" << region_end << "." << std::endl;
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

    std::cout << "Read information for " << n_samples << " samples from .fam file." << std::endl;

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

void GenoData::run_mean_imputation() {

    bool needs_impute = false;
    for (int j = 0; j < this->genotype_matrix.cols(); ++j) {
        auto col = this->genotype_matrix.col(j);
        
        if ((col.array() == -1).any()) {
            
            needs_impute = true;            

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

    if (needs_impute) {
        std::cout << "Mean imputation performed successfully." << std::endl;
    }
}

void GenoData::compute_maf() {
    maf.clear();
    mac.clear();
    maf.reserve(n_snps);
    mac.reserve(n_snps);
    for (size_t i = 0; i < n_snps; ++i) {
        const auto col = genotype_matrix.col(i);
        double alt_count = col.sum();
        double af = alt_count / (2 * n_samples);
        double maf_value = std::min(af, 1.0 - af);
        double mac_value = std::min(alt_count, 2.0 * static_cast<double>(n_samples) - alt_count);
        maf.push_back(maf_value);
        mac.push_back(static_cast<int>(std::round(mac_value)));
    }
    std::cout << "MAF computed." << std::endl;
}

void GenoData::compute_maf_problems() {
    int n_zero_maf = 0;
    int n_low_maf = 0;
    for (size_t i = 0; i < maf.size(); ++i) {
        if (std::abs(maf[i]) <= 1e-8) {
            n_zero_maf++;
        }
        if (maf[i] <= 0.05) {
            n_low_maf++;
        }
    }

    if (n_zero_maf != 0) {
        std::cout << "Warning: " << n_zero_maf << " variants with zero MAF detected." << std::endl;
    }
    if (n_low_maf != 0 && n_samples < 1000) {
        std::cout << "Warning: " << n_low_maf << " variants with MAF < 0.05 detected." << std::endl;
        std::cout << "Sample subsetting may have reduced MAF." << std::endl;
    }
}

void GenoData::filter_maf(double min_maf) {

    if (min_maf <= 0.0) {
        return;
    }

    std::vector<size_t> keep;
    keep.reserve(n_snps);
    for (size_t i = 0; i < n_snps; ++i) {
        if (maf[i] >= min_maf) {
            keep.push_back(i);
        }
    }
    if (keep.empty()) {
        std::cerr << "Error: no variants with MAF >= " << min_maf << "." << std::endl;
        exit(1);
    }
    const size_t n_drop = n_snps - keep.size();
    if (n_drop == 0) {
        std::cout << "All " << format_with_commas(n_snps) << " SNPs have MAF >= " << min_maf << "." << std::endl;
        return;
    }

    auto subset = [&keep](auto& v) {
        std::remove_reference_t<decltype(v)> out;
        out.reserve(keep.size());
        for (size_t i : keep) {
            out.push_back(v[i]);
        }
        v = std::move(out);
    };
    subset(chrom);
    subset(pos);
    subset(id);
    subset(alt);
    subset(ref);
    subset(snp_id);
    subset(index);
    subset(maf);
    subset(mac);

    Eigen::MatrixXd sliced(genotype_matrix.rows(), keep.size());
    for (size_t k = 0; k < keep.size(); ++k) {
        sliced.col(k) = genotype_matrix.col(keep[k]);
    }
    genotype_matrix = std::move(sliced);
    n_snps = keep.size();

    std::cout << "Filtered to " << format_with_commas(n_snps) << " SNPs with MAF >= " << min_maf
              << " (dropped " << format_with_commas(n_drop) << ")." << std::endl;
}

void GenoData::slice_samples(std::vector<std::string>& sample_ids) {
    
    Eigen::VectorXi rows;
    rows.resize(sample_ids.size());
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            rows(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in genotype data." << std::endl;
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

void GenoData::permute_samples(std::mt19937_64& rng) {
    const Eigen::Index n = genotype_matrix.rows();
    std::vector<Eigen::Index> source(static_cast<size_t>(n));
    std::iota(source.begin(), source.end(), static_cast<Eigen::Index>(0));
    std::shuffle(source.begin(), source.end(), rng);

    std::vector<char> placed(static_cast<size_t>(n), 0);
    for (Eigen::Index start = 0; start < n; ++start) {
        if (placed[static_cast<size_t>(start)]) {
            continue;
        }
        placed[static_cast<size_t>(start)] = 1;

        Eigen::Index here = start;
        Eigen::Index from = source[static_cast<size_t>(start)];
        while (from != start) {
            genotype_matrix.row(here).swap(genotype_matrix.row(from));
            placed[static_cast<size_t>(from)] = 1;
            here = from;
            from = source[static_cast<size_t>(from)];
        }
    }
}
