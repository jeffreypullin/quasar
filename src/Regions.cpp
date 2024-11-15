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

#include "Data.hpp"
#include "Regions.hpp"
#include "Utils.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

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

void Regions::construct_regions(FeatData& feat_data, GenoData& geno_data, int& window_size) {

    std::vector<int>& chr_f = feat_data.chrom;
    std::vector<int>& start_f = feat_data.start; 
    std::vector<int>& end_f = feat_data.end; 

    std::vector<int>& chr_g = geno_data.chrom;
    std::vector<int>& pos_g = geno_data.pos;

    int n_f = chr_f.size();
    int n_g = chr_g.size();

    std::vector<int> geno_s_ind(n_g, -1);
    std::vector<int> geno_e_ind(n_g, -1);

    std::vector<int> feat_s_ind(n_f, -1);
    std::vector<int> feat_e_ind(n_f, -1);

    int block_ind = 0;
    int i_f = 0;
    int i_f_e = 0;
    int i_g = 0;
    int i_g_e = 0;

    double magic_fraction = 0.8;

    while (i_f < n_f) {

        i_f_e = i_f;

        seek_to(i_f_e, chr_f, start_f, chr_f[i_f], end_f[i_f] + magic_fraction * window_size);
        seek_to(i_g, chr_g, pos_g, chr_f[i_f], start_f[i_f] - window_size);

        if (chr_g[i_g] == chr_f[i_f]) {

            i_g_e = i_g_e > i_g ? i_g_e : i_g;
            seek_to(i_g_e, chr_g, pos_g, chr_f[i_f], end_f[i_f_e] + window_size + 1);

            if (i_f_e >= i_f && pos_g[i_g] - end_f[i_f] <= window_size && chr_g[i_g] == chr_f[i_f]) {

                std::string region = std::to_string(chr_g[i_g]) + ":" + std::to_string(pos_g[i_g]) + "-" + std::to_string(pos_g[i_g_e]);
                regions_id.push_back(region);

                feat_s.push_back(i_f);
                feat_e.push_back(i_f_e);
                geno_s.push_back(i_g);
                geno_e.push_back(i_g_e);
                
                for (int j = i_g; j <= i_g_e; ++j){
                    if (geno_s_ind[j] == -1) {
                        geno_s_ind[j] = block_ind;
                    }
                    geno_e_ind[j] = block_ind;
                }

                for (int i = i_f; i <= i_f_e; ++i) {
                    if (feat_s_ind[i] == -1) {
                        feat_s_ind[i] = block_ind;
                    }
                    feat_e_ind[i] = block_ind;
                }
                block_ind++;
            }
        }
        i_f = i_f_e + 1;
    }
}
