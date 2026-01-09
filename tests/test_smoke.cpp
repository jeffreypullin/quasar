// Smoke tests for quasar using example data

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>
#include <cmath>

#ifndef EXAMPLE_DATA_DIR
#define EXAMPLE_DATA_DIR "."
#endif

#ifndef QUASAR_EXECUTABLE
#define QUASAR_EXECUTABLE "quasar"
#endif

int run_quasar(const std::string& args) {
    std::string cmd = std::string(QUASAR_EXECUTABLE) + " " + args + " 2>&1";
    return std::system(cmd.c_str());
}

bool file_exists(const std::string& path) {
    std::ifstream f(path);
    return f.good();
}

int count_lines(const std::string& path) {
    std::ifstream f(path);
    int count = 0;
    std::string line;
    while (std::getline(f, line)) {
        count++;
    }
    return count;
}

std::vector<std::string> split_string(const std::string& s, char delim) {
    std::vector<std::string> tokens;
    std::istringstream iss(s);
    std::string token;
    while (std::getline(iss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}

bool all_variant_pvalues_valid(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    std::getline(f, line); // skip header

    while (std::getline(f, line)) {
        auto tokens = split_string(line, '\t');
        if (tokens.size() > 9) {
            double pval = std::stod(tokens[9]);
            if (pval < 0.0 || pval > 1.0 || std::isnan(pval)) {
                return false;
            }
        }
    }
    return true;
}

bool all_gene_pvalues_valid(const std::string& path) {
    std::ifstream f(path);
    std::string line;
    std::getline(f, line); // skip header

    while (std::getline(f, line)) {
        auto tokens = split_string(line, '\t');
        if (!tokens.empty()) {
            double pval = std::stod(tokens.back());
            if (pval < 0.0 || pval > 1.0 || std::isnan(pval)) {
                return false;
            }
        }
    }
    return true;
}

void cleanup_output(const std::string& prefix) {
    std::remove((prefix + "-quasar-cis-variant.txt").c_str());
    std::remove((prefix + "-quasar-cis-region.txt").c_str());
}

TEST_CASE("LM smoke test") {
    std::string data_dir = EXAMPLE_DATA_DIR;
    std::string build_dir = BUILD_DIR;
    std::string out_prefix = build_dir + "/test_lm_output";

    cleanup_output(out_prefix);

    std::string cmd =
        "-p " + data_dir + "/chr22-n100 "
        "-b " + data_dir + "/mean-pheno-n100.bed "
        "-c " + data_dir + "/cov-n100.tsv "
        "-o " + out_prefix + " "
        "--model lm "
        "--mode cis";

    SECTION("quasar runs successfully") {
        int exit_code = run_quasar(cmd);
        REQUIRE(exit_code == 0);
    }

    SECTION("output files are created") {
        run_quasar(cmd);
        REQUIRE(file_exists(out_prefix + "-quasar-cis-variant.txt"));
        REQUIRE(file_exists(out_prefix + "-quasar-cis-region.txt"));
    }

    SECTION("variant output has expected header") {
        run_quasar(cmd);
        std::ifstream f(out_prefix + "-quasar-cis-variant.txt");
        std::string header;
        std::getline(f, header);
        REQUIRE(header.find("feature_id") != std::string::npos);
        REQUIRE(header.find("pvalue") != std::string::npos);
    }

    SECTION("all variant p-values are valid") {
        run_quasar(cmd);
        REQUIRE(all_variant_pvalues_valid(out_prefix + "-quasar-cis-variant.txt"));
    }

    SECTION("gene output has 20 genes") {
        run_quasar(cmd);
        // 1 header + 20 genes = 21 lines
        REQUIRE(count_lines(out_prefix + "-quasar-cis-region.txt") == 21);
    }

    SECTION("all gene-level p-values are valid") {
        run_quasar(cmd);
        REQUIRE(all_gene_pvalues_valid(out_prefix + "-quasar-cis-region.txt"));
    }

    cleanup_output(out_prefix);
}

TEST_CASE("NB-GLM smoke test") {
    std::string data_dir = EXAMPLE_DATA_DIR;
    std::string build_dir = BUILD_DIR;
    std::string out_prefix = build_dir + "/test_nb_glm_output";

    cleanup_output(out_prefix);

    std::string cmd =
        "-p " + data_dir + "/chr22-n100 "
        "-b " + data_dir + "/sum-pheno-n100.bed "
        "-c " + data_dir + "/cov-n100.tsv "
        "-o " + out_prefix + " "
        "--model nb_glm "
        "--use-apl "
        "--mode cis";

    SECTION("quasar runs successfully") {
        int exit_code = run_quasar(cmd);
        REQUIRE(exit_code == 0);
    }

    SECTION("output files are created") {
        run_quasar(cmd);
        REQUIRE(file_exists(out_prefix + "-quasar-cis-variant.txt"));
        REQUIRE(file_exists(out_prefix + "-quasar-cis-region.txt"));
    }

    SECTION("variant output includes phi column") {
        run_quasar(cmd);
        std::ifstream f(out_prefix + "-quasar-cis-variant.txt");
        std::string header;
        std::getline(f, header);
        REQUIRE(header.find("phi") != std::string::npos);
    }

    SECTION("all variant p-values are valid") {
        run_quasar(cmd);
        REQUIRE(all_variant_pvalues_valid(out_prefix + "-quasar-cis-variant.txt"));
    }

    SECTION("gene output has 20 genes") {
        run_quasar(cmd);
        REQUIRE(count_lines(out_prefix + "-quasar-cis-region.txt") == 21);
    }

    SECTION("all gene-level p-values are valid") {
        run_quasar(cmd);
        REQUIRE(all_gene_pvalues_valid(out_prefix + "-quasar-cis-region.txt"));
    }

    cleanup_output(out_prefix);
}

TEST_CASE("LMM smoke test") {
    std::string data_dir = EXAMPLE_DATA_DIR;
    std::string build_dir = BUILD_DIR;
    std::string out_prefix = build_dir + "/test_lmm_output";

    cleanup_output(out_prefix);

    std::string cmd =
        "-p " + data_dir + "/chr22-n100 "
        "-b " + data_dir + "/mean-pheno-n100.bed "
        "-c " + data_dir + "/cov-n100.tsv "
        "-g " + data_dir + "/grm-n100.tsv "
        "-o " + out_prefix + " "
        "--model lmm "
        "--mode cis";

    SECTION("quasar runs successfully") {
        int exit_code = run_quasar(cmd);
        REQUIRE(exit_code == 0);
    }

    SECTION("output files are created") {
        run_quasar(cmd);
        REQUIRE(file_exists(out_prefix + "-quasar-cis-variant.txt"));
        REQUIRE(file_exists(out_prefix + "-quasar-cis-region.txt"));
    }

    SECTION("all variant p-values are valid") {
        run_quasar(cmd);
        REQUIRE(all_variant_pvalues_valid(out_prefix + "-quasar-cis-variant.txt"));
    }

    SECTION("gene output has 20 genes") {
        run_quasar(cmd);
        REQUIRE(count_lines(out_prefix + "-quasar-cis-region.txt") == 21);
    }

    SECTION("all gene-level p-values are valid") {
        run_quasar(cmd);
        REQUIRE(all_gene_pvalues_valid(out_prefix + "-quasar-cis-region.txt"));
    }

    cleanup_output(out_prefix);
}
