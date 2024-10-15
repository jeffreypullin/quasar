
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
        return;
    }
    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, "\t ");
        if (tokens[0] != "feature_id") {
            std::cerr << "Error: Invalid header in phenotype file. Expected 'feature_id' as the first column." << std::endl;
            return;
        }
        sample_ids = std::vector<std::string>(tokens.begin() + 1, tokens.end());
        n_samples = sample_ids.size();
    }

    int n_pheno = 0;
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

    int row = 0;
    while (std::getline(file, line) && row < n_pheno) {
        remove_carriage_return(line);
        tokens = string_split(line, "\t ");
        if (tokens.size() != n_samples + 1) {
            std::cerr << "Error: Inconsistent number of columns in phenotype file." << std::endl;
            return;
        }
        pheno_ids.push_back(tokens[0]);
        for (int col = 0; col < n_samples; ++col) {
            data(row, col) = std::stod(tokens[col + 1]);
        }
        row++;
    }

    file.close();

    std::cout << "Read " << n_pheno << " features for " << n_samples << " samples from phenotype file." << std::endl;
}

void PhenoData::slice_samples(std::vector<std::string>& sample_ids) {
    
    Eigen::VectorXi columns;
    columns.resize(sample_ids.size());
    for (int i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            columns(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in phenotype data." << std::endl;
            return;
        }
    }
    this->data = data(Eigen::all, columns);
    this->sample_ids = sample_ids;
    n_samples = sample_ids.size();
}

void CovData::read_cov_data() {
    std::ifstream file(cov_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open covariate file: " << cov_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, "\t ");
        if (tokens.size() < 2 || tokens[0] != "sample_id") {
            std::cerr << "Error: Invalid header in covariate file. Expected 'sample_id' as the first column." << std::endl;
            return;
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

    int row = 0;
    while (std::getline(file, line) && row < n_samples) {
        tokens = string_split(line, "\t ");
        if (tokens.size() != n_cov + 1) {
            std::cerr << "Error: Inconsistent number of columns in covariate file at line " << row + 2 << std::endl;
            return;
        }
        sample_ids.push_back(tokens[0]);
        for (int col = 0; col < n_cov; ++col) {
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
    for (int i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            rows(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in phenotype data." << std::endl;
            return;
        }
    }
    this->data = data(rows, Eigen::all);
    this->sample_ids = sample_ids;
    n_samples = sample_ids.size();
}

void FeatData::read_feat_data() {

    std::ifstream file(feat_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open feature annotation file: " << feat_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, "\t ");
        if (tokens.size() < 4 || tokens[0] != "feature_id" || tokens[1] != "chrom" || tokens[2] != "start" || tokens[3] != "end") {
            std::cerr << "Error: Invalid header in feature annotation file. Expected 'feature_id', 'chrom', 'start', 'end' as the first four columns." << std::endl;
            return;
        }
    }

    n_feat = 0;
    while (std::getline(file, line)) {
        remove_carriage_return(line);
        tokens = string_split(line, "\t ");
        if (tokens.size() < 4) {
            std::cerr << "Error: Inconsistent number of columns in feature annotation file at line " << n_feat + 2 << std::endl;
            return;
        }

        feat_id.push_back(tokens[0]);
        chrom.push_back(std::stoi(tokens[1]));
        start.push_back(std::stoi(tokens[2]));
        end.push_back(std::stoi(tokens[3]));

        // Check if gene_name column exists
        if (tokens.size() > 4) {
            gene_name.push_back(tokens[4]);
        } else {
            gene_name.push_back("");
        }

        n_feat++;
    }

    file.close();

    std::cout << "Read " << n_feat << " features from feature annotation file." << std::endl;
}

void GRM::read_grm() {
    std::ifstream file(grm_file);
    if (!file.is_open()) {
        std::cerr << "Error: Unable to open GRM file: " << grm_file << std::endl;
        return;
    }

    std::string line;
    std::vector<std::string> tokens;

    if (std::getline(file, line)) {
        tokens = string_split(line, "\t ");
        if (tokens.size() < 2 || tokens[0] != "sample_id") {
            std::cerr << "Error: Invalid header in GRM file. Expected 'sample_id' as the first column." << std::endl;
            return;
        }
        n_samps = tokens.size() - 1;
        sample_ids = std::vector<std::string>(tokens.begin() + 1, tokens.end());
    }

    mat = Eigen::MatrixXd::Zero(n_samps, n_samps);

    int row = 0;
    while (std::getline(file, line)) {
        tokens = string_split(line, "\t ");
        if (tokens.size() != n_samps + 1) {
            std::cerr << "Error: Inconsistent number of columns in GRM file at line " << row + 2 << std::endl;
            return;
        }

        for (int col = 0; col < n_samps; ++col) {
            mat(row, col) = std::stod(tokens[col + 1]);
        }

        row++;
    }

    if (row != n_samps) {
        std::cerr << "Error: Number of rows does not match number of samples in GRM file." << std::endl;
        return;
    }

    file.close();

    std::cout << "Read GRM with " << n_samps << " samples." << std::endl;
}

void GRM::slice_samples(std::vector<std::string>& sample_ids) {

    Eigen::VectorXi ind;
    ind.resize(sample_ids.size());
    for (int i = 0; i < sample_ids.size(); ++i) {
        auto it = std::find(this->sample_ids.begin(), this->sample_ids.end(), sample_ids[i]);
        if (it != this->sample_ids.end()) {
            ind(i) = std::distance(this->sample_ids.begin(), it);
        } else {
            std::cerr << "Error: Sample ID " << sample_ids[i] << " not found in GRM file." << std::endl;   
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
    n_samps = sample_ids.size();
}
