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

#include "Regions.hpp"
#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void construct_regions(Param* params, Regions* regions, FeatData* feat_data, std::vector<SNP>* snps_info) {

	std::vector<int>& chr_f = feat_data->chrom;
	std::vector<int>& start_f = feat_data->start; 
	std::vector<int>& end_f = feat_data->end; 

	std::vector<int> chr_g;
	std::transform(snps_info->begin(), snps_info->end(), std::back_inserter(chr_g), 
                   [](const SNP& snp) { return snp.chrom; });

    std::vector<int> pos_g;
	std::transform(snps_info->begin(), snps_info->end(), std::back_inserter(pos_g), 
                   [](const SNP& snp) { return snp.chrom; });

	int window = params->window_size;
	
	int n_f = chr_f.size();
	int n_g = chr_g.size();
	
	std::vector<int> geno_id_s(n_g, -1);
	std::vector<int> geno_id_e(n_g, -1);
	
	std::vector<int> feat_id_s(n_f, -1);
	std::vector<int> feat_id_e(n_f, -1);
	
	int block_ind = 0;
	int i_f = 0;
	int i_f_e = 0;
	int i_g = 0;
	int i_g_e = 0;
	
	double magic_fraction = 0.8;
	
	while (i_f < n_f) {
		
		i_f_e = i_f;
		
		seek_to(i_f_e, chr_f, start_f, chr_f[i_f], end_f[i_f] + magic_fraction * window);
		seek_to(i_g, chr_g, pos_g, chr_f[i_f], start_f[i_f] - window);
		
		if (chr_g[i_g] == chr_f[i_f]) {
			
			i_g_e = i_g_e > i_g ? i_g_e : i_g;
			seek_to(i_g_e, chr_g, pos_g, chr_f[i_f], end_f[i_f_e] + window + 1);
			
			if (i_f_e >= i_f && pos_g[i_g] - end_f[i_f] <= window && chr_g[i_g] == chr_f[i_f] ){
				
				regions->feat_s.push_back(i_f);
				regions->feat_e.push_back(i_f_e);
				regions->geno_s.push_back(i_g);
				regions->geno_e.push_back(i_g_e);
				
				for(int j = i_g; j <= i_g_e; ++j){
					if (geno_id_s[j] == -1) {
                        geno_id_s[j] = block_ind;
					}
					geno_id_e[j] = block_ind;
				}
				
				std::string region = chr_g[i_g] + ":" + std::to_string(pos_g[i_g]) + "-" + std::to_string(pos_g[i_g_e]);
				
				regions->regions_id.push_back(region);

				for (int i = i_f; i <= i_f_e; ++i) {
					if (feat_id_s[i] == -1) {
						feat_id_s[i] = block_ind;
					}
					feat_id_e[i] = block_ind;
				}
				block_ind++;
			}
		}
		i_f = i_f_e + 1;
		std::cerr << block_ind << " blocks ...\n";
	}
}

void seek_to(int& i, std::vector<int>& chr, std::vector<int>& pos, int& target_chr, int target_pos){
	while (chr[i] != target_chr && i < chr.size() - 1 ){
		i++;
	}
	while (chr[i] == target_chr && pos[i] < target_pos && i < chr.size() - 1 ){
		if (chr[i + 1] == target_chr && pos[i + 1] < target_pos){
            i++;
		} else{
			break;
		}
	}
}
