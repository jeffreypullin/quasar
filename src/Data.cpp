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

#include "Data.hpp"
#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <numeric>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <set>
#include <algorithm>
#include <cmath>

void PhenoData::read_pheno_data(const std::string& mode) {

    std::ifstream file(pheno_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open phenotype file " << pheno_file << std::endl;
        exit(1);
    }
    std::string line;
    std::vector<std::string> tokens;

    size_t n_info_cols = 4;
    size_t id_col = 3;

    if (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, ",\t ");

        const bool header_has_coords =
            tokens.size() >= 4 &&
            tokens[0] == "#chr" && tokens[1] == "start" &&
            tokens[2] == "end" && tokens[3] == "phenotype_id";
        const bool header_id_only =
            tokens.size() >= 1 && tokens[0] == "phenotype_id";

        if (mode == "gwas") {
            if (header_has_coords) {
                has_genomic_coords = true;
            } else if (header_id_only) {
                has_genomic_coords = false;
                n_info_cols = 1;
                id_col = 0;
            } else {
                std::cerr << "Error: Invalid header in phenotype file. In mode 'gwas', expected either "
                          << "'#chr', 'start', 'end', 'phenotype_id' or 'phenotype_id' as the leading column names."
                          << std::endl;
                exit(1);
            }
        } else {
            has_genomic_coords = true;
            if (!header_has_coords) {
                std::cerr << "Error: Invalid header in phenotype file. Expected '#chr', 'start', 'end' and 'phenotype_id' as the first column names." << std::endl;
                exit(1);
            }
        }

        sample_ids = std::vector<std::string>(tokens.begin() + n_info_cols, tokens.end());
        n_samples = sample_ids.size();
    }

    n_pheno = 0;
    while (std::getline(file, line)) {
        if (!line.empty()) {
          n_pheno++;
        }
    }

    file.clear();
    file.seekg(0);
    // Skip header.
    std::getline(file, line);
    
    data = Eigen::MatrixXd(n_pheno, n_samples);
    pheno_ids.reserve(n_pheno);
    if (has_genomic_coords) {
        chrom.reserve(n_pheno);
        start.reserve(n_pheno);
        end.reserve(n_pheno);
    }

    size_t row = 0;
    while (std::getline(file, line) && row < n_pheno) {
        remove_carriage_return(line);
        tokens = string_split(line, ",\t ");
        if (tokens.size() != n_samples + n_info_cols) {
            std::cerr << "Error: Inconsistent number of columns in phenotype file." << std::endl;
            exit(1);
        }
        if (has_genomic_coords) {
            try {
                chrom.push_back(std::stoi(tokens[0]));
                start.push_back(std::stoi(tokens[1]));
                end.push_back(std::stoi(tokens[2]));
            } catch (const std::exception& e) {
                std::cerr << "Error: Failed to parse coordinates at row " << row + 1 
                          << " (chrom='" << tokens[0] << "', start='" << tokens[1] 
                          << "', end='" << tokens[2] << "')" << std::endl;
                exit(1);
            }
        }
        pheno_ids.push_back(tokens[id_col]);
        for (size_t col = 0; col < n_samples; ++col) {
            data(row, col) = std::stod(tokens[col + n_info_cols]);
        }
        row++;
    }

    data.transposeInPlace();

    file.close();

    std::cout << "Read " << format_with_commas(n_pheno) << " features for " 
              << format_with_commas(n_samples) << " samples from phenotype file.";
    if (mode == "gwas") {
        std::cout << " (genomic coordinates " << (has_genomic_coords ? "present" : "not provided") << ")";
    }
    std::cout << std::endl;
}

void PhenoData::prepare_sc_pheno_data() {

    std::ifstream file(pheno_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open phenotype file " << pheno_file << std::endl;
        exit(1);
    }
    std::string line;
    std::vector<std::string> tokens;

    std::getline(file, line);
    remove_carriage_return(line);
    tokens = string_split(line, ",\t ");
    if (!(tokens[0] == "sample_id" && tokens[1] == "cell_id")) {
        std::cerr << "Error: Invalid header in phenotype file. Expected 'sample_id' and 'cell_id' as the first column names." << std::endl;
        exit(1);
    }

    pheno_ids = std::vector<std::string>(tokens.begin() + 2, tokens.end());
    n_pheno = pheno_ids.size();

    sample_ids.clear();
    cell_counts.clear();
    cell_ids.clear();

    file.clear();
    file.seekg(0);
    // Skip header.
    std::getline(file, line);

    std::string prev_sample;
    bool has_prev = false;
    std::unordered_set<std::string> seen_samples;
    size_t current_count = 0;
    n_cells = 0;

    while (std::getline(file, line)) {
        remove_carriage_return(line);
        if (line.empty()) {
            continue;
        }
        tokens = string_split(line, ",\t ");
        if (tokens.size() != n_pheno + 2) {
            std::cerr << "Error: Inconsistent number of columns in phenotype file at line "
                      << n_cells + 2 << std::endl;
            exit(1);
        }

        const std::string& sample_id = tokens[0];
        if (has_prev && sample_id != prev_sample && seen_samples.count(sample_id) > 0) {
            std::cerr << "Error: Sample IDs are not contiguous in phenotype file." << std::endl;
            exit(1);
        }

        if (!has_prev) {
            prev_sample = sample_id;
            has_prev = true;
            current_count = 0;
            seen_samples.insert(sample_id);
        } else if (sample_id != prev_sample) {
            sample_ids.push_back(prev_sample);
            cell_counts.push_back(static_cast<int>(current_count));
            prev_sample = sample_id;
            current_count = 0;
            seen_samples.insert(sample_id);
        }

        current_count++;
        n_cells++;
        cell_ids.push_back(tokens[1]);
    }

    if (has_prev) {
        sample_ids.push_back(prev_sample);
        cell_counts.push_back(static_cast<int>(current_count));
    }

    n_samples = sample_ids.size();

    file.close();

    std::cout << "Detected valid single-cell data with " 
              << format_with_commas(n_pheno) << " phenotypes for "
              << format_with_commas(n_samples) << " samples and "
              << format_with_commas(n_cells) << " cells." << std::endl;
}

void PhenoData::read_anno_data(std::string anno_file) {

    std::ifstream file(anno_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open annotation file " << anno_file << std::endl;
        exit(1);
    }
    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, ",\t ");
        if (tokens.size() < 4 || tokens[0] != "#chr" || tokens[1] != "start" || tokens[2] != "end" || tokens[3] != "phenotype_id") {
            std::cerr << "Error: Invalid header in annotation file. Expected '#chr', 'start', 'end' and 'phenotype_id' as the first column names." << std::endl;
            exit(1);
        }
    }

    size_t n_pheno_file = 0;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            n_pheno_file++;
        }
    }

    std::vector<int> anno_chrom;
    std::vector<int> anno_start;
    std::vector<int> anno_end;
    std::vector<std::string> anno_pheno_ids;
    
    anno_chrom.reserve(n_pheno_file);
    anno_start.reserve(n_pheno_file);
    anno_end.reserve(n_pheno_file);
    anno_pheno_ids.reserve(n_pheno_file);

    file.clear();
    file.seekg(0);
    // Skip header.
    std::getline(file, line);

    size_t row = 0;
    while (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, ",\t ");
        try {
            anno_chrom.push_back(std::stoi(tokens[0]));
            anno_start.push_back(std::stoi(tokens[1]));
            anno_end.push_back(std::stoi(tokens[2]));
        } catch (const std::exception& e) {
            std::cerr << "Error: Failed to parse coordinates at row " << row + 1
                      << " (chrom='" << tokens[0] << "', start='" << tokens[1]
                      << "', end='" << tokens[2] << "')" << std::endl;
            exit(1);
        }
        anno_pheno_ids.push_back(tokens[3]);
        row++;
    }

    for (const auto& id : pheno_ids) {
        if (std::find(anno_pheno_ids.begin(), anno_pheno_ids.end(), id) == anno_pheno_ids.end()) {
            std::cerr << "Error: phenotype_id '" << id << "' not found in annotation file." << std::endl;
        }
    }
    // Filter chrom, start, end to only the entries in pheno_ids, using anno_pheno_ids as guide
    std::unordered_map<std::string, size_t> anno_idx_map;
    for (size_t i = 0; i < anno_pheno_ids.size(); ++i) {
        anno_idx_map[anno_pheno_ids[i]] = i;
    }

    chrom.resize(pheno_ids.size());
    start.resize(pheno_ids.size());
    end.resize(pheno_ids.size());

    for (size_t i = 0; i < pheno_ids.size(); ++i) {
        const auto& pid = pheno_ids[i];
        auto it = anno_idx_map.find(pid);
        size_t idx = it->second;
        chrom[i] = anno_chrom[idx];
        start[i] = anno_start[idx];
        end[i] = anno_end[idx];
    }

    file.close();

    std::cout << "Read " << format_with_commas(n_pheno_file) << " features from annotation file, " << 
        "retaining information for " << format_with_commas(n_pheno) << " features." << std::endl;
}

void PhenoData::filter_pheno_ids(int filt_chrom) {
    std::vector<std::string> new_pheno_ids;
    std::vector<int> new_chrom;
    std::vector<int> new_start;
    std::vector<int> new_end;

    for (size_t i = 0; i < n_pheno; ++i) {
        if (chrom[i] == filt_chrom) {
            new_pheno_ids.push_back(pheno_ids[i]);
            new_chrom.push_back(chrom[i]);
            new_start.push_back(start[i]);
            new_end.push_back(end[i]);
            pheno_inds.push_back(i);
        }
    }

    if (new_pheno_ids.empty()) {
        std::cerr << "Error: no phenotypes found on chromosome " << filt_chrom << "." << std::endl;
        exit(1);
    }

    pheno_ids = new_pheno_ids;
    chrom = new_chrom;
    start = new_start;
    end = new_end;
    n_pheno = pheno_ids.size();

    std::cout << "Filtered to " << format_with_commas(n_pheno) << " phenotypes on chromosome " << filt_chrom << "." << std::endl;
}

void PhenoData::read_sc_pheno_data() {

    std::ifstream file(pheno_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open phenotype file " << pheno_file << std::endl;
        exit(1);
    }
    std::string line;
    std::vector<std::string> tokens;

    std::getline(file, line);
    remove_carriage_return(line);
    tokens = string_split(line, ",\t ");
    size_t total_pheno_cols = tokens.size() - 2;

    sc_data = Eigen::MatrixXd(n_cells, n_pheno);
    offset.resize(n_cells);

    size_t row = 0;
    while (std::getline(file, line)) {
        remove_carriage_return(line);
        if (line.empty()) {
            continue;
        }
        tokens = string_split(line, ",\t ");

        double sum = 0.0;
        for (size_t i = 2; i < total_pheno_cols + 2; ++i) {
            sum += std::stod(tokens[i]);
        }
        if (sum <= 0.0) {
            offset(static_cast<Eigen::Index>(row)) = 0.0;
        } else {
            offset(static_cast<Eigen::Index>(row)) = std::log(sum);
        }

        for (size_t j = 0; j < pheno_inds.size(); ++j) {
            size_t col_idx = pheno_inds[j] + 2;
            sc_data(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(j)) =
                std::stod(tokens[col_idx]);
        }
        row++;
    }

    file.close();

    std::cout << "Read " << format_with_commas(n_pheno) << " phenotypes for "
        << format_with_commas(n_cells) << " cells in phenotype file." << std::endl;
}

void PhenoData::write_pheno_data(std::string out_file) {
    std::ofstream file(out_file);

    if (has_genomic_coords) {
        file << "#chr\tstart\tend\tphenotype_id";
    } else {
        file << "phenotype_id";
    }
    for (const auto& sample_id : sample_ids) {
        file << "\t" << sample_id;
    }
    file << "\n";

    Eigen::MatrixXd transposed_data = data.transpose();

    for (size_t i = 0; i < n_pheno; ++i) {
        if (has_genomic_coords) {
            file << chrom[i] << "\t" << start[i] << "\t" << end[i] << "\t" << pheno_ids[i];
        } else {
            file << pheno_ids[i];
        }
        for (size_t j = 0; j < n_samples; ++j) {
            file << "\t" << transposed_data(i, j);
        }
        file << "\n";
    }

    file.close();
}

void PhenoData::slice_samples(std::vector<std::string>& sample_ids) {
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
    Eigen::MatrixXd temp = this->data;
    this->data = temp(rows, Eigen::all);
    this->sample_ids = sample_ids;
    n_samples = sample_ids.size();
}

void PhenoData::slice_sc_samples(std::vector<std::string>& sample_ids) {

    if (this->sample_ids.size() != this->cell_counts.size()) {
        std::cerr << "Error: sample_ids and cell_counts have inconsistent sizes in phenotype data." << std::endl;
        exit(1);
    }

    std::unordered_map<std::string, size_t> sample_index;
    sample_index.reserve(this->sample_ids.size());
    for (size_t i = 0; i < this->sample_ids.size(); ++i) {
        sample_index[this->sample_ids[i]] = i;
    }

    std::vector<size_t> block_starts(this->sample_ids.size());
    size_t running = 0;
    for (size_t i = 0; i < this->cell_counts.size(); ++i) {
        block_starts[i] = running;
        running += static_cast<size_t>(this->cell_counts[i]);
    }

    std::vector<int> new_cell_counts;
    new_cell_counts.reserve(sample_ids.size());
    size_t new_n_cells = 0;
    for (const auto& sample_id : sample_ids) {
        auto it = sample_index.find(sample_id);
        if (it == sample_index.end()) {
            std::cerr << "Error: Sample ID " << sample_id << " not found in phenotype data." << std::endl;
            exit(1);
        }
        int count = this->cell_counts[it->second];
        new_cell_counts.push_back(count);
        new_n_cells += static_cast<size_t>(count);
    }

    Eigen::MatrixXd new_sc_data(new_n_cells, n_pheno);
    Eigen::VectorXd new_offset(new_n_cells);
    std::vector<std::string> new_cell_ids;
    new_cell_ids.reserve(new_n_cells);

    size_t out_row = 0;
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        size_t idx = sample_index[sample_ids[i]];
        size_t start = block_starts[idx];
        int count = this->cell_counts[idx];
        if (count > 0) {
            new_sc_data.block(out_row, 0, count, n_pheno) = sc_data.block(start, 0, count, n_pheno);
            new_offset.segment(out_row, count) = offset.segment(start, count);
            if (!cell_ids.empty()) {
                new_cell_ids.insert(new_cell_ids.end(),
                                    cell_ids.begin() + start,
                                    cell_ids.begin() + start + static_cast<size_t>(count));
            }
            out_row += static_cast<size_t>(count);
        }
    }

    this->sc_data = new_sc_data;
    this->offset = new_offset;
    this->sample_ids = sample_ids;
    this->cell_counts = new_cell_counts;
    if (!cell_ids.empty()) {
        this->cell_ids = std::move(new_cell_ids);
    } else {
        this->cell_ids.clear();
    }
    n_samples = sample_ids.size();
    n_cells = new_n_cells;
}

void PhenoData::slice_chromosome(int chrom_id) {

    if (chrom.size() == 0) {
        std::cout << "Error: chromosome not initialised." << std::endl;
        exit(1);
    }
    // Assume window parameters are not initialised.
    std::vector<int> new_chrom;
    std::vector<int> new_start;
    std::vector<int> new_end;
    std::vector<std::string> new_pheno_ids;
    std::vector<int> col_inds;
    for (size_t i = 0; i < n_pheno; ++i) {
        if (chrom[i] == chrom_id) {
            new_chrom.push_back(chrom[i]);
            new_start.push_back(start[i]);
            new_end.push_back(end[i]);
            new_pheno_ids.push_back(pheno_ids[i]); 
            col_inds.push_back(i); 
        }
    }

    if (data_type != "single-cell") {
        Eigen::VectorXi cols;
        cols.resize(col_inds.size());
        for (size_t i = 0; i < col_inds.size(); ++i) {
            cols(i) = col_inds[i];
        }
        Eigen::MatrixXd new_data(data.rows(), col_inds.size());
        for (size_t i = 0; i < col_inds.size(); ++i) {
            if (col_inds[i] >= data.cols()) {
                std::cerr << "Error: Invalid column index " << col_inds[i] << std::endl;
                exit(1);
            }
            new_data.col(i) = data.col(col_inds[i]);
        }
        this->data = new_data;
    }

    this->chrom = new_chrom;
    this->start = new_start;
    this->end = new_end;
    this->pheno_ids = new_pheno_ids;
    this->n_pheno = new_chrom.size();
}

void PhenoData::construct_windows(GenoData& geno_data, int window_size, bool verbose) {

    for (size_t i = 0; i < n_pheno; ++i) {

        int window_start = 0;
        int window_end = 0; 
        int window_n = 0;
        
        int chr_f = chrom[i];
        int start_pos_f = start[i];
        int end_pos_f = end[i];

        std::vector<int> g_chr_vec = geno_data.chrom;
        std::vector<int> g_pos_vec = geno_data.pos;
        for (size_t j = 0; j < geno_data.n_snps - 1; ++j) {
            if (g_chr_vec[j] == chr_f) {

                if (window_start == 0 && g_pos_vec[j] >= start_pos_f - window_size) {
                    window_start = j;
                }
                if (g_pos_vec[j] > end_pos_f + window_size) {
                    window_end = j;
                    break;
                }
            }
        }

        if (window_start == 0 && window_end == 0) {
            if (verbose && g_chr_vec[0] == chr_f) {
                std::cout << "Warning: No variants found in window for feature " << pheno_ids[i] << std::endl;
            }
        } else if (window_end == 0) {
            int index_last_on_chr = -1;
            for (size_t k = 0; k < geno_data.n_snps - 1; ++k) {
                if (g_chr_vec[k] == chr_f) {
                    index_last_on_chr = k;
                }
            }
            window_end = index_last_on_chr;
        }

        window_n = window_end - window_start;

        this->window_start.push_back(window_start);
        this->window_end.push_back(window_end);
        this->window_n.push_back(window_n);
    }
}

void CovData::read_cov_data() {
    std::ifstream file(cov_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open covariate file: " << cov_file << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, ",\t ");
        if (tokens.size() < 2 || tokens[0] != "sample_id") {
            std::cerr << "Error: Invalid header in covariate file. Expected 'sample_id' as the first column." << std::endl;
            exit(1);
        }
        cov_ids = std::vector<std::string>(tokens.begin() + 1, tokens.end());
        n_cov = cov_ids.size();
    }

    n_samples = 0;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            n_samples++;
        }
    }

    file.clear();
    file.seekg(0);
    // Skip header.
    std::getline(file, line);

    data = Eigen::MatrixXd(n_samples, n_cov);
    sample_ids.reserve(n_samples);

    size_t row = 0;
    while (std::getline(file, line) && row < n_samples) {
        tokens = string_split(line, ",\t ");
        if (tokens.size() != n_cov + 1) {
            std::cerr << "Error: Inconsistent number of columns in covariate file at line " << row + 2 << std::endl;
            exit(1);
        }
        sample_ids.push_back(tokens[0]);
        for (size_t col = 0; col < n_cov; ++col) {
            data(row, col) = std::stod(tokens[col + 1]);
        }
        row++;
    }

    file.close();

    // Add intercept column if not already present in the data.
    bool has_intercept = false;
    for (size_t j = 0; j < n_cov; j++) {
        bool is_intercept = true;
        double first_val = data(0, j);
        for (size_t i = 1; i < n_samples; i++) {
            if (data(i, j) != first_val) {
                is_intercept = false;
                break;
            }
        }
        if (is_intercept && first_val == 1.0) {
            has_intercept = true;
            break;
        }
    }

    if (!has_intercept) {
        Eigen::MatrixXd new_data(n_samples, n_cov + 1);
        new_data << Eigen::VectorXd::Ones(n_samples), data;
        data = new_data;
        cov_ids.insert(cov_ids.begin(), "intercept");
        n_cov++;
    }

    std::cout << "Read " << n_cov << " covariates for "<< format_with_commas(n_samples) << " samples from covariate file." << std::endl;
}

void CovData::check_cov_data_type() {
    std::ifstream file(cov_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open covariate file: " << cov_file << std::endl;
        exit(1);
    }

    std::string line;
    if (std::getline(file, line)) {
        std::vector<std::string> tokens = string_split(line, ",\t ");
        if (tokens.size() > 1 && tokens[0] == "sample_id" && tokens[1] == "cell_id") {
            cov_data_type = "single-cell";
        } else if (tokens.size() > 0 && tokens[0] == "sample_id") {
            cov_data_type = "bulk";
        } else {
            std::cerr << "Error: Covariate file header must start with 'sample_id' or with 'sample_id, cell_id'" << std::endl;
            exit(1);
        }
    }
    file.close();

}

bool CovData::is_covariate_categorical() {
    static const size_t max_unique_values = 10;
    
    const Eigen::MatrixXd* active_cov_data = &data;
    if (cov_data_type == "single-cell" && sc_data.size() > 0) {
        active_cov_data = &sc_data;
    }

    std::set<double> unique_values;
    for (Eigen::Index i = 0; i < active_cov_data->rows(); ++i) {
        double x = (*active_cov_data)(i, interaction_ind);
        unique_values.insert(x);
        if (unique_values.size() > max_unique_values) {
            return false;
        }
    }
    return true;
}

void CovData::add_squared_covariate() {
    std::string sq_covariate_id = interaction_id + "_sq";

    if (cov_data_type == "single-cell" && sc_data.size() > 0) {
        Eigen::MatrixXd updated_sc_data(sc_data.rows(), sc_data.cols() + 1);
        updated_sc_data.leftCols(sc_data.cols()) = sc_data;
        updated_sc_data.col(sc_data.cols()) = sc_data.col(interaction_ind).array().square().matrix();
        sc_data = updated_sc_data;
    } else {
        Eigen::MatrixXd updated_data(data.rows(), data.cols() + 1);
        updated_data.leftCols(data.cols()) = data;
        updated_data.col(data.cols()) = data.col(interaction_ind).array().square().matrix();
        data = updated_data;
    }

    cov_ids.push_back(sq_covariate_id);
    n_cov++;
}

void CovData::add_bw_covariates() {
    Eigen::VectorXd int_col = sc_data.col(interaction_ind);
    int n_cells_total = sc_data.rows();
    int old_ncols = static_cast<int>(sc_data.cols());

    Eigen::VectorXd b_col = Eigen::VectorXd::Zero(n_cells_total);
    int row = 0;
    for (size_t i = 0; i < cell_counts.size(); ++i) {
        int nc = cell_counts[i];
        double mean_i = int_col.segment(row, nc).mean();
        b_col.segment(row, nc).setConstant(mean_i);
        row += nc;
    }
    Eigen::VectorXd w_col = int_col - b_col;

    // Drop the original interaction column and append the between/within replacements.
    Eigen::MatrixXd updated(n_cells_total, old_ncols - 1 + 2);
    int dest = 0;
    for (int j = 0; j < old_ncols; ++j) {
        if (j == interaction_ind) continue;
        updated.col(dest++) = sc_data.col(j);
    }
    updated.col(dest)     = b_col;
    updated.col(dest + 1) = w_col;
    sc_data = updated;

    cov_ids.erase(cov_ids.begin() + interaction_ind);
    cov_ids.push_back(interaction_id + "_b");
    cov_ids.push_back(interaction_id + "_w");
    n_cov += 1;

    interaction_id  = interaction_id + "_w";
    interaction_ind = static_cast<int>(sc_data.cols()) - 1;
}

void CovData::standardisze_data() {
    Eigen::MatrixXd* matrix = &data;
    if (cov_data_type == "single-cell") {
        matrix = (sc_data.size() > 0) ? &sc_data : &data;
    } else if (sc_data.size() > 0) {
        matrix = &sc_data;
    }

    if (matrix->size() == 0) {
        return;
    }

    for (Eigen::Index col = 0; col < matrix->cols(); ++col) {
        double mean = matrix->col(col).mean();
        Eigen::ArrayXd centered = matrix->col(col).array() - mean;
        double var = centered.square().mean();
        if (var <= 1e-8) {
            continue;
        }
        matrix->col(col).array() = centered / std::sqrt(var);
    }
}

void CovData::read_sc_cov_data() {
    std::ifstream file(cov_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open covariate file: " << cov_file << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, ",\t ");
        if (tokens.size() < 3 || tokens[0] != "sample_id" || tokens[1] != "cell_id") {
            std::cerr << "Error: Invalid header in covariate file. Expected 'sample_id' and 'cell_id' as the first columns." << std::endl;
            exit(1);
        }
        cov_ids = std::vector<std::string>(tokens.begin() + 2, tokens.end());
        n_cov = cov_ids.size();
    }

    n_cells = 0;
    while (std::getline(file, line)) {
        if (!line.empty()) {
            n_cells++;
        }
    }

    file.clear();
    file.seekg(0);
    // Skip header.
    std::getline(file, line);

    data = Eigen::MatrixXd(n_cells, n_cov);
    sample_ids.clear();
    cell_ids.clear();
    cell_counts.clear();
    sample_ids.reserve(n_cells);
    cell_ids.reserve(n_cells);

    std::string prev_sample;
    bool has_prev = false;
    std::unordered_set<std::string> seen_samples;
    size_t current_count = 0;

    size_t row = 0;
    while (std::getline(file, line) && row < n_cells) {
        remove_carriage_return(line);
        if (line.empty()) {
            continue;
        }
        tokens = string_split(line, ",\t ");
        const std::string& sample_id = tokens[0];
        if (has_prev && sample_id != prev_sample && seen_samples.count(sample_id) > 0) {
            std::cerr << "Error: Sample IDs are not contiguous in covariate file." << std::endl;
            exit(1);
        }
        if (!has_prev) {
            prev_sample = sample_id;
            has_prev = true;
            seen_samples.insert(sample_id);
            current_count = 0;
        } else if (sample_id != prev_sample) {
            // Finalize previous sample block.
            sample_ids.push_back(prev_sample);
            cell_counts.push_back(static_cast<int>(current_count));
            prev_sample = sample_id;
            seen_samples.insert(sample_id);
            current_count = 0;

        }

        cell_ids.push_back(tokens[1]);
        current_count++;
        for (size_t col = 0; col < n_cov; ++col) {
            data(row, col) = std::stod(tokens[col + 2]);
        }
        row++;
    }

    if (has_prev) {
        sample_ids.push_back(prev_sample);
        cell_counts.push_back(static_cast<int>(current_count));
    }

    n_samples = sample_ids.size();

    file.close();

    // Add intercept column if not already present in the data.
    bool has_intercept = false;
    for (size_t j = 0; j < n_cov; j++) {
        bool is_intercept = true;
        double first_val = data(0, j);
        for (size_t i = 1; i < n_cells; i++) {
            if (data(i, j) != first_val) {
                is_intercept = false;
                break;
            }
        }
        if (is_intercept && first_val == 1.0) {
            has_intercept = true;
            break;
        }
    }

    if (!has_intercept) {
        Eigen::MatrixXd new_data(n_cells, n_cov + 1);
        new_data << Eigen::VectorXd::Ones(n_cells), data;
        sc_data = new_data;
        cov_ids.insert(cov_ids.begin(), "intercept");
        n_cov++;
    }

    std::cout << "Read " << n_cov << " covariates for "
              << format_with_commas(n_cells) << " cells ("
              << format_with_commas(n_samples) << " samples) from covariate file." << std::endl;
}

void CovData::expand_cov_data(std::vector<int> cell_counts) {
    
    n_cells = std::accumulate(cell_counts.begin(), cell_counts.end(), 0);

    Eigen::MatrixXd expanded_data;
    expanded_data.resize(n_cells, n_cov);

    int row_idx = 0;
    for (size_t i = 0; i < cell_counts.size(); ++i) {
        for (int j = 0; j < cell_counts[i]; ++j) {
            expanded_data.row(row_idx) = data.row(i);
            row_idx++;
        }
    }

    sc_data = expanded_data;
}

void CovData::slice_samples(std::vector<std::string>& sample_ids) {
    
    Eigen::VectorXi rows;
    rows.resize(sample_ids.size());
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            rows(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in covariate data." << std::endl;
            exit(1);
        }
    }
    Eigen::MatrixXd temp = this->data;
    this->data = temp(rows, Eigen::all);
    this->sample_ids = sample_ids;
    n_samples = sample_ids.size();
}

void CovData::slice_sc_samples(std::vector<std::string>& sample_ids) {

    std::unordered_map<std::string, size_t> sample_index;
    sample_index.reserve(this->sample_ids.size());
    for (size_t i = 0; i < this->sample_ids.size(); ++i) {
        sample_index[this->sample_ids[i]] = i;
    }

    // Compute block starts from cell_counts (prefix sums).
    std::vector<size_t> block_starts(this->cell_counts.size());
    size_t running = 0;
    for (size_t i = 0; i < this->cell_counts.size(); ++i) {
        block_starts[i] = running;
        running += static_cast<size_t>(this->cell_counts[i]);
    }

    std::vector<int> new_cell_counts;
    new_cell_counts.reserve(sample_ids.size());
    size_t new_n_cells = 0;
    for (const auto& sid : sample_ids) {
        auto it = sample_index.find(sid);
        if (it == sample_index.end()) {
            std::cerr << "Error: Sample ID " << sid << " not found in covariate data." << std::endl;
            exit(1);
        }
        int count = this->cell_counts[it->second];
        new_cell_counts.push_back(count);
        new_n_cells += static_cast<size_t>(count);
    }

    Eigen::MatrixXd new_sc_data(new_n_cells, n_cov);
    std::vector<std::string> new_cell_ids;
    new_cell_ids.reserve(new_n_cells);

    size_t out_row = 0;
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        size_t idx = sample_index[sample_ids[i]];
        size_t bstart = block_starts[idx];
        int count = this->cell_counts[idx];
        if (count > 0) {
            new_sc_data.block(out_row, 0, count, n_cov) = sc_data.block(bstart, 0, count, n_cov);
            new_cell_ids.insert(new_cell_ids.end(),
                                cell_ids.begin() + bstart,
                                cell_ids.begin() + bstart + static_cast<size_t>(count));
            out_row += static_cast<size_t>(count);
        }
    }

    this->sc_data = new_sc_data;
    this->sample_ids = sample_ids;
    this->cell_counts = new_cell_counts;
    this->cell_ids = std::move(new_cell_ids);

    n_samples = sample_ids.size();
    n_cells = new_n_cells;
}

size_t align_sc_cell_ids(PhenoData& pheno_data, CovData& cov_data) {

    const size_t pheno_before = pheno_data.n_cells;

    if (pheno_data.cell_ids.size() != pheno_data.n_cells || cov_data.cell_ids.size() != cov_data.n_cells) {
        std::cerr << "Error: cell_ids length does not match n_cells; cannot align cell_ids." << std::endl;
        exit(1);
    }

    // Compute block starts (prefix sums).
    std::vector<size_t> pheno_starts(pheno_data.cell_counts.size());
    std::vector<size_t> cov_starts(cov_data.cell_counts.size());
    size_t run_p = 0, run_c = 0;
    for (size_t i = 0; i < pheno_data.cell_counts.size(); ++i) {
        pheno_starts[i] = run_p;
        cov_starts[i] = run_c;
        run_p += static_cast<size_t>(pheno_data.cell_counts[i]);
        run_c += static_cast<size_t>(cov_data.cell_counts[i]);
    }

    std::vector<int> new_cell_counts;
    new_cell_counts.reserve(pheno_data.sample_ids.size());
    size_t kept_total = 0;

    for (size_t i = 0; i < pheno_data.sample_ids.size(); ++i) {
        const size_t p_start = pheno_starts[i];
        const int p_count = pheno_data.cell_counts[i];
        const size_t c_start = cov_starts[i];
        const int c_count = cov_data.cell_counts[i];

        std::unordered_map<std::string, int> cov_pos;
        // Factor of 2 to keep the hash tables load factor low.
        cov_pos.reserve(static_cast<size_t>(c_count) * 2 + 1);

        for (int j = 0; j < c_count; ++j) {
            const std::string& cid = cov_data.cell_ids[c_start + static_cast<size_t>(j)];
            auto ins = cov_pos.emplace(cid, j);
            if (!ins.second) {
                std::cerr << "Error: Duplicate cell_id '" << cid
                          << "' within sample '" << cov_data.sample_ids[i]
                          << "' in covariate data." << std::endl;
                exit(1);
            }
        }

        int kept = 0;
        for (int j = 0; j < p_count; ++j) {
            const std::string& cid = pheno_data.cell_ids[p_start + static_cast<size_t>(j)];
            if (cov_pos.find(cid) != cov_pos.end()) {
                kept++;
            }
        }
        new_cell_counts.push_back(kept);
        kept_total += static_cast<size_t>(kept);
    }

    if (kept_total == 0) {
        std::cerr << "Error: No overlapping cell_ids between phenotype and covariate data after alignment." << std::endl;
        exit(1);
    }

    // Second pass: construct aligned/filtered matrices in phenotype order.
    Eigen::MatrixXd new_pheno_sc(kept_total, pheno_data.n_pheno);
    Eigen::VectorXd new_offset(kept_total);
    Eigen::MatrixXd new_cov(kept_total, cov_data.n_cov);
    std::vector<std::string> new_cell_ids;
    new_cell_ids.reserve(kept_total);

    size_t out = 0;
    for (size_t i = 0; i < pheno_data.sample_ids.size(); ++i) {
        const size_t p_start = pheno_starts[i];
        const int p_count = pheno_data.cell_counts[i];
        const size_t c_start = cov_starts[i];
        const int c_count = cov_data.cell_counts[i];

        std::unordered_map<std::string, int> cov_pos;
        cov_pos.reserve(static_cast<size_t>(c_count) * 2 + 1);
        for (int j = 0; j < c_count; ++j) {
            cov_pos.emplace(cov_data.cell_ids[c_start + static_cast<size_t>(j)], j);
        }

        for (int j = 0; j < p_count; ++j) {
            const size_t p_row = p_start + static_cast<size_t>(j);
            const std::string& cid = pheno_data.cell_ids[p_row];
            auto it = cov_pos.find(cid);
            if (it == cov_pos.end()) {
                continue;
            }
            const size_t c_row = c_start + static_cast<size_t>(it->second);

            new_pheno_sc.row(out) = pheno_data.sc_data.row(static_cast<Eigen::Index>(p_row));
            new_offset(static_cast<Eigen::Index>(out)) = pheno_data.offset(static_cast<Eigen::Index>(p_row));
            new_cov.row(out) = cov_data.sc_data.row(static_cast<Eigen::Index>(c_row));
            new_cell_ids.push_back(cid);
            out++;
        }
    }

    if (out != kept_total) {
        std::cerr << "Error: Internal alignment error (kept_total mismatch)." << std::endl;
        exit(1);
    }

    pheno_data.sc_data = new_pheno_sc;
    pheno_data.offset = new_offset;
    pheno_data.cell_ids = new_cell_ids;
    pheno_data.cell_counts = new_cell_counts;
    pheno_data.n_cells = kept_total;

    cov_data.sc_data = new_cov;
    cov_data.cell_ids = new_cell_ids;
    cov_data.cell_counts = new_cell_counts;
    cov_data.n_cells = kept_total;

    double percent_retained = (pheno_before == 0) ? 0.0 : (static_cast<double>(kept_total) / pheno_before) * 100.0;
    std::cout << "Aligned cell IDs between phenotype and covariate data.\n"
              << "Retained " << std::fixed << std::setprecision(0) << percent_retained << "%"
              << " of phenotype cells." << std::endl;
    std::cout << std::defaultfloat << std::setprecision(6);

    return kept_total;
}

void CellGroups::read_cell_groups() {

    std::ifstream f(file);
    if (!f.is_open()) {
        std::cerr << "Error: Unable to open cell groups file: " << file << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string> tokens;

    if (!std::getline(f, line)) {
        std::cerr << "Error: cell groups file is empty: " << file << std::endl;
        exit(1);
    }
    remove_carriage_return(line);
    tokens = string_split(line, ",\t ");
    if (tokens.size() < 2 || tokens[0] != "group" || tokens[1] != "cell_id") {
        std::cerr << "Error: Invalid header in cell groups file. Expected 'group' and 'cell_id' as the first two columns." << std::endl;
        exit(1);
    }
    has_values = (tokens.size() >= 3 && tokens[2] == "value");

    std::unordered_map<std::string, int> group_id_to_idx;
    cell_id_to_group_idx_.clear();
    cell_id_to_value_.clear();
    group_ids.clear();

    size_t row = 1;
    while (std::getline(f, line)) {
        row++;
        remove_carriage_return(line);
        if (line.empty()) {
            continue;
        }
        tokens = string_split(line, ",\t ");
        if (tokens.size() < 2) {
            std::cerr << "Error: Cell groups file at line " << row << " has fewer than 2 columns." << std::endl;
            exit(1);
        }
        const std::string& group = tokens[0];
        const std::string& cell_id = tokens[1];

        auto inserted = cell_id_to_group_idx_.emplace(cell_id, -1);
        if (!inserted.second) {
            std::cerr << "Error: Duplicate cell_id '" << cell_id
                      << "' in cell groups file at line " << row << "." << std::endl;
            exit(1);
        }

        if (has_values) {
            if (tokens.size() < 3) {
                std::cerr << "Error: Cell groups file at line " << row << " is missing the 'value' column." << std::endl;
                exit(1);
            }
            try {
                cell_id_to_value_[cell_id] = std::stod(tokens[2]);
            } catch (...) {
                std::cerr << "Error: Cell groups file at line " << row << " has a non-numeric value '" << tokens[2] << "'." << std::endl;
                exit(1);
            }
        }

        auto git = group_id_to_idx.find(group);
        int gi;
        if (git == group_id_to_idx.end()) {
            gi = static_cast<int>(group_ids.size());
            group_id_to_idx.emplace(group, gi);
            group_ids.push_back(group);
        } else {
            gi = git->second;
        }
        inserted.first->second = gi;
    }

    n_groups = group_ids.size();
    f.close();

    if (n_groups == 0) {
        std::cerr << "Error: No groups found in cell groups file." << std::endl;
        exit(1);
    }

    std::cout << "Read " << format_with_commas(cell_id_to_group_idx_.size())
              << " cell-to-group assignments across " << n_groups
              << " groups from cell groups file." << std::endl;
}

void CellGroups::align_to_cells(const std::vector<std::string>& cell_ids) {

    if (n_groups == 0) {
        std::cerr << "Error: align_to_cells called before read_cell_groups." << std::endl;
        exit(1);
    }

    cell_to_group.assign(cell_ids.size(), -1);
    cells_per_group.assign(n_groups, std::vector<size_t>());

    size_t n_assigned = 0;
    for (size_t i = 0; i < cell_ids.size(); ++i) {
        auto it = cell_id_to_group_idx_.find(cell_ids[i]);
        if (it == cell_id_to_group_idx_.end()) {
            std::cerr << "Error: cell_id '" << cell_ids[i]
                      << "' is present in phenotype data but missing from cell groups file." << std::endl;
            exit(1);
        }
        int gi = it->second;
        cell_to_group[i] = gi;
        cells_per_group[gi].push_back(i);
        n_assigned++;
    }

    if (n_assigned != cell_id_to_group_idx_.size()) {
        std::cerr << "Warning: " << (cell_id_to_group_idx_.size() - n_assigned)
                  << " cell_id(s) in cell groups file are not present in phenotype data; ignoring."
                  << std::endl;
    }

    for (size_t gi = 0; gi < n_groups; ++gi) {
        if (cells_per_group[gi].empty()) {
            std::cerr << "Error: cell group '" << group_ids[gi]
                      << "' has no cells in the aligned phenotype data." << std::endl;
            exit(1);
        }
    }

    if (has_values) {
        group_values.resize(n_groups);
        for (size_t gi = 0; gi < n_groups; ++gi) {
            std::vector<double> vals;
            vals.reserve(cells_per_group[gi].size());
            for (size_t cell_idx : cells_per_group[gi]) {
                const std::string& cid = cell_ids[cell_idx];
                auto vit = cell_id_to_value_.find(cid);
                if (vit != cell_id_to_value_.end()) {
                    vals.push_back(vit->second);
                }
            }
            if (vals.empty()) {
                group_values[gi] = std::numeric_limits<double>::quiet_NaN();
            } else {
                std::sort(vals.begin(), vals.end());
                size_t mid = vals.size() / 2;
                group_values[gi] = (vals.size() % 2 == 0)
                    ? (vals[mid - 1] + vals[mid]) / 2.0
                    : vals[mid];
            }
        }
        cell_id_to_value_.clear();
    }
}

void GRM::read_grm() {
    std::ifstream file(grm_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open GRM file: " << grm_file << std::endl;
        exit(1);
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, ",\t ");
        if (tokens.size() < 2 || tokens[0] != "sample_id") {
            std::cerr << "Error: Invalid header in GRM file. Expected 'sample_id' as the first column." << std::endl;
            exit(1);
        }
        n_samples = tokens.size() - 1;
        sample_ids = std::vector<std::string>(tokens.begin() + 1, tokens.end());
    }

    mat = Eigen::MatrixXd::Zero(n_samples, n_samples);

    size_t row = 0;
    while (std::getline(file, line)) {
        tokens = string_split(line, ",\t ");
        if (tokens.size() != n_samples + 1) {
            std::cerr << "Error: Inconsistent number of columns in GRM file at line " << row + 2 << std::endl;
            exit(1);
        }

        for (size_t col = 0; col < n_samples; ++col) {
            mat(row, col) = std::stod(tokens[col + 1]);
        }

        row++;
    }

    if (row != n_samples) {
        std::cerr << "Error: Number of rows does not match number of samples in GRM file." << std::endl;
        exit(1);
    }

    file.close();

    std::cout << "Read GRM with " << format_with_commas(n_samples) << " samples." << std::endl;
}

void GRM::slice_samples(std::vector<std::string>& sample_ids) {

    Eigen::VectorXi ind;
    ind.resize(sample_ids.size());
    for (size_t i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            ind(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in GRM file." << std::endl;
            exit(1); 
            return;
        }
    }
    Eigen::MatrixXd sliced_mat(ind.size(), ind.size());
    for (int i = 0; i < ind.size(); ++i) {
        for (int j = 0; j < ind.size(); ++j) {
            sliced_mat(i, j) = mat(ind(i), ind(j));
        }
    }
    this->mat = sliced_mat;

    this->sample_ids = sample_ids;
    n_samples = sample_ids.size();
}
