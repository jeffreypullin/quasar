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

#ifndef DATA_H
#define DATA_H

#include <unordered_map>

#include "Quasar.hpp"

class PhenoData {

    public:
      std::string pheno_file;
      std::string data_type;
      
      size_t n_pheno;
      size_t n_samples;
      size_t n_cells;

      std::vector<std::string> pheno_ids;
      std::vector<size_t> pheno_inds;
      std::vector<std::string> sample_ids;
      std::vector<std::string> cell_ids;
      std::vector<int> cell_counts;

      Eigen::MatrixXd data;
      Eigen::MatrixXd sc_data;
      Eigen::VectorXd offset;

      std::vector<int> chrom;
      std::vector<int> start;
      std::vector<int> end;
      std::vector<int> window_start;
      std::vector<int> window_end;
      std::vector<int> window_n;

      bool has_genomic_coords = true;

      PhenoData(std::string pheno_file, std::string data_type) {
        this->pheno_file = pheno_file;
        this->data_type = data_type;
      }

      void read_pheno_data(const std::string& mode);
      void write_pheno_data(std::string out_file);
      
      void construct_windows(GenoData& geno_data, int window_size, bool verbose); 

      void slice_chromosome(int chrom_id);
      void slice_samples(std::vector<std::string>& sample_ids);

      // Single-cell specific.
      void prepare_sc_pheno_data();
      void filter_pheno_ids(int filt_chrom);
      void read_sc_pheno_data(bool compute_offset);
      void slice_sc_samples(std::vector<std::string>& sample_ids);
      void read_anno_data(std::string anno_file);
};

class CovData {

    public:
      std::string cov_file;
      std::string data_type;
      std::string cov_data_type;

      size_t n_cov;
      size_t n_samples;
      size_t n_cells;

      std::vector<std::string> cov_ids;
      std::vector<std::string> sample_ids;
      std::vector<std::string> cell_ids;
      std::vector<int> cell_counts;
      Eigen::MatrixXd data;
      Eigen::MatrixXd sc_data;

      std::vector<int> interaction_inds;
      std::vector<std::string> interaction_ids;

      CovData(std::string cov_file) {
        this->cov_file = cov_file;
      }
    
      void check_cov_data_type();
      void read_cov_data();
      void read_sc_cov_data();
      void expand_cov_data(std::vector<int> cell_counts);
      void collapse_cov_data();
      void add_bw_covariates();
      void standardisze_data();

      void slice_samples(std::vector<std::string>& sample_ids);
      void slice_sc_samples(std::vector<std::string>& sample_ids);
};

class GRM {

  public:
    std::string grm_file;
    size_t n_samples;
    std::vector<std::string> sample_ids;
    Eigen::MatrixXd mat;
    GRM(std::string grm_file) {
      this->grm_file = grm_file;
    }
    void slice_samples(std::vector<std::string>& sample_ids);
    void read_grm();
};

size_t align_sc_cell_ids(PhenoData& pheno_data, CovData& cov_data);

class CellGroups {

  public:
    std::string file;
    size_t n_groups = 0;
    std::vector<std::string> group_ids;
    std::vector<int> cell_to_group;
    std::vector<std::vector<size_t>> cells_per_group;
    bool has_values = false;
    std::vector<double> group_values;

    CellGroups(std::string file) {
      this->file = file;
    }

    void read_cell_groups();
    void align_to_cells(const std::vector<std::string>& cell_ids);

  private:
    std::unordered_map<std::string, int> cell_id_to_group_idx_;
    std::unordered_map<std::string, double> cell_id_to_value_;
};

class OffsetData {

  public:
    std::string file;
    std::string offset_data_type;

    OffsetData(std::string file) {
      this->file = file;
    }

    void read_offset_data(const std::string& data_type);
    Eigen::VectorXd align_to_samples(const std::vector<std::string>& sample_ids);
    Eigen::VectorXd align_to_cells(const std::vector<std::string>& cell_ids);

  private:
    std::unordered_map<std::string, double> key_to_offset_;
};

#endif