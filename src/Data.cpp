
#include "Data.hpp"
#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void PhenoData::read_pheno_data() {

    std::ifstream file(pheno_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open phenotype file " << pheno_file << std::endl;
        exit(1);
    }
    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, ",\t ");
        if (tokens[0] != "#chr" || tokens[1] != "start" || 
            tokens[2] != "end" || tokens[3] != "phenotype_id") {
            std::cerr << "Error: Invalid header in phenotype file. Expected '#chr', 'start', 'end' and 'phenotype_id' as the first column names." << std::endl;
            exit(1);
        }
        sample_ids = std::vector<std::string>(tokens.begin() + 4, tokens.end());
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
    chrom.reserve(n_pheno);
    start.reserve(n_pheno);
    end.reserve(n_pheno);

    size_t row = 0;
    while (std::getline(file, line) && row < n_pheno) {
        remove_carriage_return(line);
        tokens = string_split(line, ",\t ");
        if (tokens.size() != n_samples + 4) {
            std::cerr << "Error: Inconsistent number of columns in phenotype file." << std::endl;
            exit(1);
        }
        chrom.push_back(std::stoi(tokens[0]));
        start.push_back(std::stoi(tokens[1]));
        end.push_back(std::stoi(tokens[2]));
        std::cout << "Start: " << tokens[1] << ", End: " << tokens[2] << std::endl;
        pheno_ids.push_back(tokens[3]);
        for (size_t col = 0; col < n_samples; ++col) {
            data(row, col) = std::stod(tokens[col + 4]);
        }
        row++;
    }

    data.transposeInPlace();

    file.close();

    std::cout << "Read " << n_pheno << " features for " << n_samples << " samples from phenotype file." << std::endl;
}

void PhenoData::write_pheno_data(std::string out_file) {
    std::ofstream file(out_file);
    
    file << "#chr\tstart\tend\tphenotype_id";
    for (const auto& sample_id : sample_ids) {
        file << "\t" << sample_id;
    }
    file << "\n";

    Eigen::MatrixXd transposed_data = data.transpose();

    for (size_t i = 0; i < n_pheno; ++i) {
        file << chrom[i] << "\t" << start[i] << "\t" << end[i] << "\t" << pheno_ids[i];
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

    Eigen::VectorXi cols;
    cols.resize(col_inds.size());
    for (size_t i = 0; i < col_inds.size(); ++i) {
        cols(i) = col_inds[i];
    }
    Eigen::MatrixXd new_data(data.rows(), col_inds.size());
    for (size_t i = 0; i < col_inds.size(); ++i) {
        if (col_inds[i] >= data.cols()) {
            std::cerr << "Error: Invalid row index " << col_inds[i] << std::endl;
            exit(1);
        }
        new_data.col(i) = data.col(col_inds[i]);
    }
    this->data = new_data;

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
            window_end = geno_data.n_snps - 1;
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

    std::cout << "Read " << n_cov << " covariates for "<< n_samples << " samples from covariate file." << std::endl;
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

    std::cout << "Read GRM with " << n_samples << " samples." << std::endl;
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
