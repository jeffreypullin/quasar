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

#ifndef REGIONS_HPP
#define REGIONS_HPP

#include <vector>
#include <string>

#include "Data.hpp"
#include "Quasar.hpp"
#include "Geno.hpp"

struct Regions {    

    // Here s denotes start and e denotes end.
    std::vector<int> feat_s;
    std::vector<int> feat_e;
    std::vector<int> geno_s;
    std::vector<int> geno_e;
	std::vector<int> geno_s_ind;
	std::vector<int> geno_e_ind;
	std::vector<int> feat_s_ind;
	std::vector<int> feat_e_ind;
	std::vector<std::string> regions_id;

};

void construct_regions(Param* params, Regions* regions, FeatData* feat_data, std::vector<SNP>* snps_info);
void seek_to(int& i, std::vector<int>& chr, std::vector<int>& pos, int& target_chr, int target_pos);

#endif