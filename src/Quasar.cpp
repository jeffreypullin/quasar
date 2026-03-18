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

#include "Quasar.hpp"
#include "Data.hpp"
#include "ModelFit.hpp"
#include "Residualise.hpp"
#include "ScoreTest.hpp"
#include "Geno.hpp"
#include "Utils.hpp"
#include <iostream>
#include <utility>

int main(int argc, char* argv[]) {

    Params params;

    cxxopts::Options options("quasar", "QTL mapping software");
    options.add_options()
        ("h,help", "Display help message")
        ("v,version", "Display version information")
        // Data arguments.
        ("p,plink", "Prefix to PLINK files (.bed, .bim, .fam)", cxxopts::value<std::string>(params.plink_prefix))
        ("b,bed", "Bed file holding phenotype informaton", cxxopts::value<std::string>(params.bed_file)->default_value("no-bed"))
        ("sc-pheno", "File holding single-cell level phenotype data", cxxopts::value<std::string>(params.sc_pheno_file))
        ("anno", "File holding feature annotation", cxxopts::value<std::string>(params.anno_file))
        ("c,cov", "Covariate file", cxxopts::value<std::string>(params.cov_file))
        ("r,resid", "Residualised phenotype bed file", cxxopts::value<std::string>(params.resid_file)->default_value("no-resid"))
        ("f,fit", "Model fit file", cxxopts::value<std::string>(params.fit_file)->default_value("no-fit"))
        ("g,grm", "Genomic relatedness matrix", cxxopts::value<std::string>(params.grm_file)->default_value("no-grm"))
        ("i,interaction", "Covariate column name for GxE interaction testing", cxxopts::value<std::string>(params.interaction_cov))
        // Execution arguments.
        ("mode", "Mode to run quasar in (residualise, cis, trans, gwas)", cxxopts::value<std::string>(params.mode))
        ("model", "Statistical model to use for QTL mapping (lmm, glmm)", cxxopts::value<std::string>(params.model))
        ("w,window", "Cis window size in base pairs", cxxopts::value<int>(params.window_size))
        ("use-apl", "Use adjusted profile likelihood to estimate NB dispersion", cxxopts::value<bool>(params.use_apl))
        ("use-quant-res", "Use randomised quantile residuals to compute score tests with NB-GLM model", cxxopts::value<bool>(params.use_quant_res))
        ("pheno-chr", "Only map QTLs for genes on this chromosome", cxxopts::value<int>(params.pheno_chr))
        // Output arguments.
        ("o,out", "Output file prefix", cxxopts::value<std::string>(params.out))
        ("verbose", "Run with extensive output to terminal", cxxopts::value<bool>(params.verbose));

    // Parse the arguments.
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (result.count("version")) {
        std::cout << "quasar version 1.1.0" << std::endl;
        exit(0);
    }

    if (params.out == "") {
        params.out = "quasar_output";
    }
    params.do_interaction = !params.interaction_cov.empty();

    std::cout << "\nquasar execution started." << std::endl;

    if (params.model != "lmm" && 
        params.model != "p_glmm" && 
        params.model != "lm" && 
        params.model != "p_glm" && 
        params.model != "nb_glm" && 
        params.model != "p_glmm_id" &&
        params.model != "p_glmm_sc" &&
        params.model != "p_glmm_grm" &&
        params.model != "nb_glmm") {
        std::cerr << "Invalid model specified. Please use 'lm', 'lmm', 'p_glm', 'nb_glm', 'p_glmm', "
                  << "'p_glmm_id', 'p_glmm_sc', 'p_glmm_grm' or 'nb_glmm'." << std::endl;
        exit(1);
    }

    if (params.mode != "cis" && params.mode != "trans" && params.mode != "gwas" && params.mode != "residualise") {
        std::cerr << "Invalid mode specified. Please use one of 'cis', 'trans', 'gwas', 'residualise'" << std::endl;
        exit(1);
    }

    if (params.model != "nb_glm" && (params.use_apl)) {
        std::cerr << "Error: The --use-apl flagged is only applicable when usisng the NB-GLM model." << std::endl;
        exit(1);
    }
    
    if (params.model != "nb_glm" && (params.use_quant_res)) {
        std::cerr << "Error: The --use-quant-res flagged is only applicable when usisng the NB-GLM model." << std::endl;
        exit(1);
    }

    if (result.count("pheno-chr")) {
        if (params.pheno_chr < 1 || params.pheno_chr > 22) {
            std::cerr << "Error: --pheno-chr argument must integer between 1 and 22." << std::endl;
            exit(1);
        }
    }

    if (!params.sc_pheno_file.empty()) {
        params.data_type = "single-cell";
        if (params.anno_file.empty()) {
            std::cerr << "Error: feature annotation file (--anno) must also be specified for single-cell data." << std::endl;
            exit(1);
        }
        if (params.model != "p_glmm_sc") {
            std::cerr << "Error: only the `p_glmm_sc` model is compatible with single-cell data." << std::endl;
            exit(1);
        }
    } else {
        params.data_type = "bulk";
    }
    
    std::cout << "\nMode: " << params.mode << std::endl;
    std::cout << "Model: " << params.model << std::endl;
    std::cout << "Data type: " << params.data_type << std::endl;
    if (params.do_interaction) {
        std::cout << "\nPerforming interaction testing" << std::endl;
        std::cout << "Interaction covariate: " << params.interaction_cov << std::endl;
    }

    if (params.model == "p_glm") {
        std::cout << "\nWarning: using the Poisson GLM is not recommended due to its high rate of false positives." << std::endl;
    }

    if (params.model == "nb_glmm") {
        std::cout << "\nWarning: using the NB-GLMM is not recommended, use the Poisson GLMM instead." << std::endl;
    }

    bool mixed_model = params.model == "lmm" || params.model == "p_glmm" || params.model == "nb_glmm";
    if (!mixed_model && params.grm_file != "no-grm") {
       std::cout << "\nA GRM is not needed when using the LM, P-GLM or NB-GLM models and will be ignored." << std::endl;
    }

    if (params.mode == "residualise" && params.resid_file != "no-resid") {
        std::cerr << "\nError: cannot use residuals as input in mode 'residualise'" << std::endl;
        exit(1);
    }

    std::cout << "\nReading genotype data information..." << std::endl;
    GenoData geno_data(params.plink_prefix);
    geno_data.read_fam_file();
    geno_data.read_bim_file();

    bool one_chrom = false;
    if (params.mode == "cis" || params.mode == "residualise") {
        std::vector<int> g_chrom = geno_data.chrom;
        one_chrom = std::equal(g_chrom.begin() + 1, g_chrom.end(), g_chrom.begin());
        if (one_chrom) {
            std::cout << "\nOnly one chromosome detected in genotype data." << std::endl;
        }
    }

    std::cout << "\nReading non-genotype data..." << std::endl;
    std::string pheno_file;


    if (params.data_type == "single-cell") {
        if (params.sc_pheno_file.empty()) {
            std::cerr << "Error: --sc_pheno argument must be specified when data type is single-cell." << std::endl;
            exit(1);
        } else {
            pheno_file = params.sc_pheno_file;
        }
    }

    bool use_resid = false;
    if (params.data_type == "bulk") {
        bool has_bed = !params.bed_file.empty() && params.bed_file != "no-bed";
        bool has_resid = !params.resid_file.empty() && params.resid_file != "no-resid";
        if (has_bed && has_resid) {
            std::cerr << "Error: Only one of --bed or --resid can be specified for bulk data type." << std::endl;
            exit(1);
        }
        if (!has_bed && !has_resid) {
            std::cerr << "Error: One of --bed or --resid must be specified for bulk data type." << std::endl;
            exit(1);
        }
        if (has_bed) {
            pheno_file = params.bed_file;
        } else if (has_resid) {
            pheno_file = params.resid_file;
            use_resid = true;
        }
    }

    PhenoData pheno_data(pheno_file, params.data_type);
    if (params.data_type == "single-cell") {
        pheno_data.prepare_sc_pheno_data();
        pheno_data.read_anno_data(params.anno_file);
        if (result.count("pheno-chr")) {
            pheno_data.filter_pheno_ids(params.pheno_chr);
        } else if (one_chrom && params.mode == "cis") {
            pheno_data.filter_pheno_ids(geno_data.chrom.front());
        }
        pheno_data.read_sc_pheno_data();
    } else {
        pheno_data.read_pheno_data();
    }
    
    CovData cov_data(params.cov_file);
    cov_data.check_cov_data_type();
    if (cov_data.cov_data_type == "single-cell") {
        cov_data.read_sc_cov_data();
    } else {
        cov_data.read_cov_data();
    }

    if (params.do_interaction) {
        auto it = std::find(cov_data.cov_ids.begin(), cov_data.cov_ids.end(), params.interaction_cov);
        if (it == cov_data.cov_ids.end()) {
            std::cerr << "Error: interaction covariate '" << params.interaction_cov
                      << "' not found in covariate file columns." << std::endl;
            exit(1);
        }
    }

    GRM grm(params.grm_file);
    if (mixed_model) {
        if (params.grm_file == "no-grm") {
            std::cerr << "Error: GRM file must be provided when using mixed models." << std::endl;
            exit(1);
        }
        grm.read_grm();
    }

    std::cout << "\nComputing sample intersection and filtering data..." << std::endl;
    std::vector<std::vector<std::string>> sample_ids_vecs = {
        geno_data.sample_ids,
        pheno_data.sample_ids, 
        cov_data.sample_ids, 
    };
    if (mixed_model) {
        sample_ids_vecs.push_back(grm.sample_ids);
    }
    std::vector<std::string> int_sample_ids = intersection(sample_ids_vecs);

    if (int_sample_ids.size() == 0) {
        std::cerr << "Error: no common sample ids found." << std::endl;
        exit(1);
    }

    if (params.data_type == "single-cell") {
        pheno_data.slice_sc_samples(int_sample_ids);
    } else {
        pheno_data.slice_samples(int_sample_ids);
    }

    if (cov_data.cov_data_type == "single-cell") {
        cov_data.slice_sc_samples(int_sample_ids);
    } else {
        cov_data.slice_samples(int_sample_ids);
    }

    if (params.data_type == "single-cell" && cov_data.cov_data_type == "single-cell") {
        align_sc_cell_ids(pheno_data, cov_data);
    }
    if (mixed_model) {
        grm.slice_samples(int_sample_ids);
    }
    std::cout << "Running analysis for " << int_sample_ids.size() << " common samples across data inputs." << std::endl;

    if (params.do_interaction) {
        cov_data.interaction_id = params.interaction_cov;
        cov_data.interaction_ind = std::distance(
            cov_data.cov_ids.begin(),
            std::find(
                cov_data.cov_ids.begin(),
                cov_data.cov_ids.end(),
                params.interaction_cov
            )
        );
    }

    if (params.data_type == "bulk") {
        std::vector<int> g_chrom = geno_data.chrom;
        bool one_chrom = std::equal(g_chrom.begin() + 1, g_chrom.end(), g_chrom.begin());
        if (result.count("pheno-chr")) {
            std::cout << "Filtering phenotype data to features on chromosome: " << params.pheno_chr << std::endl;
            pheno_data.slice_chromosome(params.pheno_chr);
        } else if (one_chrom) {
            std::cout << "\nMode 'cis' and only one chromosome detected." << std::endl;
            std::cout << "Filtering phenotype data to features on chromosome: " << g_chrom.front() << std::endl;
            pheno_data.slice_chromosome(g_chrom.front());
        }
    }
    std::cout << "\nRunning analysis for " << format_with_commas(pheno_data.n_pheno) << " phenotypes." << std::endl;

    ModelFit model_fit(params.model, params.fit_file, pheno_data);
    if (!use_resid) {
        
        // Check model and data align.
        // FIXME: Add check for single cell data.
        if (params.model == "p_glmm" ||
            params.model == "p_glm" ||
            params.model == "nb_glm" || 
            params.model == "nb_glmm" ||
            params.model == "p_glmm_grm"
        ) {     
            Eigen::VectorXd first_gene = pheno_data.data.col(0).head(10);
            bool has_negative = (first_gene.array() < 0).any();
            bool has_noninteger = ((first_gene.array() - first_gene.array().floor()) > 0).any();
        
            if (has_negative || has_noninteger) {
                std::cerr << "Error: GLM/GLMM models require count data. The first gene appears to contain non-count values." << std::endl;
                exit(1);
            }
        }

        std::cout << "\nResidualising data..." << std::endl;
        residualise(params, model_fit, cov_data, pheno_data, grm);
        std::cout << "\nResidualisation finished." << std::endl;
    } else {
        model_fit.read_model_fit();
    }

    if (params.mode == "residualise") {
        std::cout << "\nWriting residuals to file..." << std::endl;
        pheno_data.write_pheno_data(params.out + "-resids.bed");
        std::cout << "Residuals written to file." << std::endl;

        std::cout << "\nWriting model fit to file..." << std::endl;
        model_fit.write_model_fit(params.out);
        std::cout << "Model fit written to file." << std::endl;

        std::cout << "\nquasar execution finished." << std::endl;
        exit(0);
    }

    std::cout << "\nReading genotype data..." << std::endl;
    geno_data.prepare_bed_file();
    geno_data.read_bed_file();
    geno_data.run_mean_imputation();
    geno_data.slice_samples(int_sample_ids);
    geno_data.compute_maf();
    geno_data.compute_maf_problems();

    std::cout << "\nConstructing cis-windows..." << std::endl;
    pheno_data.construct_windows(geno_data, params.window_size, params.verbose);
    std::cout << "Cis-windows constructed." << std::endl;

    std::cout << "\nPerforming variant score tests..." << std::endl;
    score_test(params, model_fit, geno_data, pheno_data, cov_data);
    std::cout << "Variant score tests finished." << std::endl;

    std::cout << "\nquasar execution finished." << std::endl;
    return 0;
}
