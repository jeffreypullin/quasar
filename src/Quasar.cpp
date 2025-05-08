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

#include "Quasar.hpp"
#include "Data.hpp"
#include "Geno.hpp"
#include "QTLMapping.hpp"
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
        ("b,bed", "Bed file holding phenotype informaton", cxxopts::value<std::string>(params.bed_file))
        ("c,cov", "Covariate file", cxxopts::value<std::string>(params.cov_file))
        ("g,grm", "Genomic relatedness matrix", cxxopts::value<std::string>(params.grm_file)->default_value("no-grm"))
        // Model arguments.
        ("m,model", "Statistical model to use for QTL mapping (lmm, glmm)", cxxopts::value<std::string>(params.model))
        ("w,window", "Cis window size in base pairs", cxxopts::value<int>(params.window_size))
        ("use-apl", "Use adjusted profile likelihood to estimate NB dispersion", cxxopts::value<bool>(params.use_apl))
        // Output arguments.
        ("o,output-prefix", "Output file prefix", cxxopts::value<std::string>(params.output_prefix))
        ("verbose", "Run with extensive output to terminal", cxxopts::value<bool>(params.verbose));

    // Parse the arguments.
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (result.count("version")) {
        std::cout << "quasar version 0.1.0" << std::endl;
        exit(0);
    }

    if (params.output_prefix == "") {
        params.output_prefix = "quasar_output";
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

    bool mixed_model = params.model == "lmm" || params.model == "p_glmm" || params.model == "nb_glmm";
    if (!mixed_model && params.grm_file != "no-grm") {
       std::cout << "\nA GRM is not needed when using the LM or NB-GLM models and will be ignored." << std::endl;
    }

    std::cout << "\nReading genotype data..." << std::endl;
    GenoData geno_data(params.plink_prefix);
    geno_data.read_bim_file();
    geno_data.read_fam_file();
    geno_data.prepare_bed_file();
    geno_data.read_bed_file();
    geno_data.run_mean_imputation();

    std::cout << "\nReading non-genotype data..." << std::endl;
    PhenoData pheno_data(params.bed_file);
    pheno_data.read_pheno_data();
    CovData cov_data(params.cov_file);
    cov_data.read_cov_data();
    GRM grm(params.grm_file);
    if (mixed_model) {
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

    geno_data.slice_samples(int_sample_ids);
    pheno_data.slice_samples(int_sample_ids);
    cov_data.slice_samples(int_sample_ids);
    if (mixed_model) {
        grm.slice_samples(int_sample_ids);
    }

    std::cout << "Running analysis for " << int_sample_ids.size() << " common samples across data inputs" << std::endl;

    // Check model and data align.
    if (params.model == "glm" || params.model == "glmm") {
        Eigen::VectorXd first_gene = pheno_data.data.col(0).head(10);
        bool has_negative = (first_gene.array() < 0).any();
        bool has_noninteger = ((first_gene.array() - first_gene.array().floor()) > 0).any();
        
        if (has_negative || has_noninteger) {
            std::cerr << "Error: GLM/GLMM models require count data. The first gene appears to contain non-count values." << std::endl;
            exit(1);
        }
    }

    std::cout << "\nFiltering phenotype data..." << std::endl;

    std::vector<int> g_chrom = geno_data.chrom;
    bool one_chrom = std::equal(g_chrom.begin() + 1, g_chrom.end(), g_chrom.begin());
    if (one_chrom) {
        std::cout << "Only one chromsome detected." << std::endl;
        std::cout << "Filtering phenotype data to chromosome: " << g_chrom.front() << std::endl;
        pheno_data.slice_chromosome(g_chrom.front());
    }
    std::cout << "\nRunning analysis for " << pheno_data.n_pheno << " phenotypes." << std::endl;
    
    std::cout << "\nConstructing cis-windows..." << std::endl;
    pheno_data.construct_windows(geno_data, params.window_size, params.verbose);
    std::cout << "Cis-windows constructed." << std::endl;

    std::cout << "\nRunning QTL mapping..." << std::endl;
    run_qtl_mapping(params, geno_data, cov_data, pheno_data, grm);
    std::cout << "\nQTL mapping finished." << std::endl;

    std::cout << "\nquasar execution finished." << std::endl;
    return 0;
}
