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
        ("c,cov", "Covariate file", cxxopts::value<std::string>(params.cov_file))
        ("r,resid", "Residualised phenotype bed file", cxxopts::value<std::string>(params.resid_file)->default_value("no-resid"))
        ("f,fit", "Model fit file", cxxopts::value<std::string>(params.fit_file)->default_value("no-fit"))
        ("g,grm", "Genomic relatedness matrix", cxxopts::value<std::string>(params.grm_file)->default_value("no-grm"))
        // Execution arguments.
        ("mode", "Mode to run quasar in (residualise, cis, trans, gwas)", cxxopts::value<std::string>(params.mode))
        ("model", "Statistical model to use for QTL mapping (lmm, glmm)", cxxopts::value<std::string>(params.model))
        ("w,window", "Cis window size in base pairs", cxxopts::value<int>(params.window_size))
        ("use-apl", "Use adjusted profile likelihood to estimate NB dispersion", cxxopts::value<bool>(params.use_apl))
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
        std::cout << "quasar version 1.0.0" << std::endl;
        exit(0);
    }

    if (params.out == "") {
        params.out = "quasar_output";
    }

    std::cout << "\nquasar execution started." << std::endl;

    // Check model.
    if (params.model != "lmm" && 
        params.model != "p_glmm" && 
        params.model != "lm" && 
        params.model != "p_glm" && 
        params.model != "nb_glm" && 
        params.model != "nb_glmm") {
        std::cerr << "Invalid model specified. Please use 'lm', 'lmm', 'p_glm', 'nb_glm', 'p_glmm', 'nb_glmm'." << std::endl;
        exit(1);
    }

    // Check mode.
    if (params.mode != "cis" && params.mode != "trans" && params.mode != "gwas" && params.mode != "residualise") {
        std::cerr << "Invalid mode specified. Please use one of 'cis', 'trans', 'gwas', 'residualise'" << std::endl;
        exit(1);
    }
    
    std::cout << "\nmode: " << params.mode << std::endl;
    std::cout << "model: " << params.model << std::endl;

    if (params.model == "p_glm") {
        std::cout << "Warning: using the Poisson GLM is not recommended due to its high rate of false positives." << std::endl;
    }

    if (params.model == "nb_glmm") {
        std::cout << "Warning: using the NB-GLMM is not recommended, use the Poisson GLMM instead." << std::endl;
    }

    bool mixed_model = params.model == "lmm" || params.model == "p_glmm" || params.model == "nb_glmm";
    if (!mixed_model && params.grm_file != "no-grm") {
       std::cout << "\nA GRM is not needed when using the LM, P-GLM or NB-GLM models and will be ignored." << std::endl;
    }

    if (params.mode == "trans" || params.mode == "gwas") {
        if (params.resid_file == "no-resid" || params.fit_file == "no-fit") {
           std::cerr << "\nError: the residual and model fit files must be specified for modes 'trans' and 'gwas'." << std::endl; 
           exit(1);
        }
    }

    if (params.mode == "residualise" && params.resid_file != "no-resid") {
        std::cerr << "\nError: cannot use residuals as input in mode 'residualise'" << std::endl;
        exit(1);
    }

    std::cout << "\nReading non-genotype data..." << std::endl;

    std::string pheno_file;
    if (params.mode == "cis" || params.mode == "residualise") {
        pheno_file = params.bed_file;
    } else {
        pheno_file = params.resid_file;
    }
    PhenoData pheno_data(pheno_file);
    pheno_data.read_pheno_data();
    
    CovData cov_data(params.cov_file);
    cov_data.read_cov_data();
    
    GRM grm(params.grm_file);
    if (mixed_model) {
        grm.read_grm();
    }
    
    GenoData geno_data(params.plink_prefix);
    geno_data.read_fam_file();
    geno_data.read_bim_file();

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

    pheno_data.slice_samples(int_sample_ids);
    cov_data.slice_samples(int_sample_ids);
    if (mixed_model) {
        grm.slice_samples(int_sample_ids);
    }
    std::cout << "Running analysis for " << int_sample_ids.size() << " common samples across data inputs." << std::endl;

    if (params.mode == "cis") {
        std::vector<int> g_chrom = geno_data.chrom;
        bool one_chrom = std::equal(g_chrom.begin() + 1, g_chrom.end(), g_chrom.begin());
        if (one_chrom) {
            std::cout << "\nMode 'cis' and only one chromosome detected." << std::endl;
            std::cout << "Filtering phenotype data to features on chromosome: " << g_chrom.front() << std::endl;
            pheno_data.slice_chromosome(g_chrom.front());
        }
    }
    std::cout << "\nRunning analysis for " << format_with_commas(pheno_data.n_pheno) << " phenotypes." << std::endl;

    ModelFit model_fit(params.model, params.fit_file, pheno_data);
    if (params.mode == "cis" || params.mode == "residualise") {
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

    std::cout << "\nConstructing cis-windows..." << std::endl;
    pheno_data.construct_windows(geno_data, params.window_size, params.verbose);
    std::cout << "Cis-windows constructed." << std::endl;

    std::cout << "\nPerforming variant score tests..." << std::endl;
    score_test(params, model_fit, geno_data, pheno_data, cov_data);
    std::cout << "Variant score tests finished." << std::endl;

    std::cout << "\nquasar execution finished." << std::endl;
    return 0;
}
