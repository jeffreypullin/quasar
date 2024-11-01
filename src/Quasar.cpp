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
#include "Regions.hpp"
#include "QTLMapping.hpp"
#include "Utils.hpp"
#include <iostream>
#include <utility>

int main(int argc, char* argv[]) {

    Param params;

    cxxopts::Options options("quasar", "QTL mapping software");
    options.add_options()
        ("h,help", "Display help message")
        ("v,version", "Display version information")
        ("m,model", "Statistical model to use for QTL mapping (lmm, glmm)", cxxopts::value<std::string>(params.model))
        ("f,feat", "Feature annotation file", cxxopts::value<std::string>(params.feat_file))
        ("c,cov", "Covariate file", cxxopts::value<std::string>(params.cov_file))
        ("p,pheno", "Phenotype file", cxxopts::value<std::string>(params.pheno_file))
        ("b,bed", "Prefix to PLINK files (.bed, .bim, .fam)", cxxopts::value<std::string>(params.bed_prefix))
        ("g,grm", "Genomic relatedness matrix", cxxopts::value<std::string>(params.grm_file))
        ("o,output-prefix", "Output directory prefix", cxxopts::value<std::string>(params.output_prefix))
        ("w,window", "Cis window size in base pairs", cxxopts::value<int>(params.window_size));

    // Parse the arguments.
    auto result = options.parse(argc, argv);

    // Handle arguments
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if (result.count("version")) {
        std::cout << "quasar version 0.0.1" << std::endl;
        exit(0);
    }

    std::cout << "\nquasar execution started." << std::endl;

    std::cout << "\nReading genotype data..." << std::endl;
    GenoData geno_data(params.bed_prefix);
    geno_data.read_bim_file();
    geno_data.read_fam_file();
    geno_data.prepare_bed_file();
    geno_data.read_bed_file();
    geno_data.compute_variant_var();

    std::cout << "\nReading non-genotype data..." << std::endl;
    PhenoData pheno_data(params.pheno_file);
    pheno_data.read_pheno_data();
    CovData cov_data(params.cov_file);
    cov_data.read_cov_data();
    FeatData feat_data(params.feat_file);
    feat_data.read_feat_data();
    GRM grm(params.grm_file);
    grm.read_grm();

    std::cout << "\nComputing sample intersection and filtering data..." << std::endl;
    std::vector<std::vector<std::string>> sample_ids_vecs = {
        grm.sample_ids, 
        cov_data.sample_ids, 
        pheno_data.sample_ids, 
        geno_data.sample_ids
    };
    std::vector<std::string> int_sample_ids = intersection(sample_ids_vecs);

    pheno_data.slice_samples(int_sample_ids);
    cov_data.slice_samples(int_sample_ids);
    grm.slice_samples(int_sample_ids);
    geno_data.slice_samples(int_sample_ids);
     
    std::cout << "Running analysis for " << int_sample_ids.size() << " common samples across data inputs" << std::endl;

    std::cout << "\nConstructing regions..." << std::endl;
    Regions regions(feat_data, geno_data, params.window_size);
    std::cout << "Regions constructed." << std::endl;
    
    std::cout << "\nRunning QTL mapping..." << std::endl;
    std::cout << "Using model: " << params.model << std::endl;
    if (params.model == "lmm") {
        run_qtl_mapping_lmm(geno_data, feat_data, cov_data, pheno_data, grm, regions);
    } else if (params.model == "glmm") {
        run_qtl_mapping_glmm(geno_data, feat_data, cov_data, pheno_data, grm, regions);
    } else {
        std::cerr << "Invalid model specified. Please use 'lmm' or 'glmm'." << std::endl;
        exit(1);
    }
    std::cout << "\nQTL mapping finished." << std::endl;

    std::cout << "\nquasar execution finished." << std::endl;
    return 0;
}
