# Apache License Version 2.0
# Copyright 2022 Xin Huang
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import allel
import math
import numpy as np

#@profile
def parse_ind_file(filename):
    """
    Description:
        Helper function to read sample information from files.

    Arguments:
        filename str: Name of the file containing sample information.

    Returns:
        samples list: Sample information.
    """
  
    f = open(filename, 'r') 
    samples = [l.rstrip() for l in f.readlines()] 
    f.close()

    if len(samples) == 0:
        raise Exception(f'No sample is found in {filename}! Please check your data.')

    return samples

#@profile
def read_geno_data(vcf, ind, anc_allele_file, filter_missing):
    """
    Description:
        Helper function to read genotype data from VCF files.

    Arguments:
        vcf str: Name of the VCF file containing genotype data.
        ind list: List containing names of samples.
        anc_allele_file str: Name of the BED file containing ancestral allele information.
        filter_missing bool: Indicating whether filtering missing data or not.

    Returns:
        data dict: Genotype data.
    """

    vcf = allel.read_vcf(vcf, alt_number=1, samples=ind)
    gt = vcf['calldata/GT']
    chr_names = np.unique(vcf['variants/CHROM'])
    samples = vcf['samples']
    pos = vcf['variants/POS']
    ref = vcf['variants/REF']
    alt = vcf['variants/ALT']

    if anc_allele_file != None: anc_allele = read_anc_allele(anc_allele_file)
    data = dict()
    for c in chr_names:
        if c not in data.keys():
            data[c] = dict()
            data[c]['POS'] = pos
            data[c]['REF'] = ref
            data[c]['ALT'] = alt
            data[c]['GT'] = gt
        index = np.where(vcf['variants/CHROM'] == c)
        data = filter_data(data, c, index)
        # Remove missing data
        if filter_missing:
            index = data[c]['GT'].count_missing(axis=1) == len(samples)
            data = filter_data(data, c, ~index)
        if anc_allele_file != None: data = check_anc_allele(data, anc_allele, c)

    return data

#@profile
def filter_data(data, c, index):
    """
    Description:
        Helper function to filter genotype data.

    Arguments:
        data dict: Genotype data for filtering.
            c str: Names of chromosomes.
        index numpy.ndarray: A boolean array determines variants to be removed.

    Returns:
        data dict: Genotype data after filtering.
    """

    data[c]['POS'] = data[c]['POS'][index]
    data[c]['REF'] = data[c]['REF'][index]
    data[c]['ALT'] = data[c]['ALT'][index]
    data[c]['GT'] = allel.GenotypeArray(data[c]['GT'][index])

    return data

#@profile
def read_data(vcf, ref_ind_file, tgt_ind_file, src_ind_file, anc_allele_file):
    """
    Description:
        Helper function for reading data from reference and target populations.

    Arguments:
        vcf str: Name of the VCF file containing genotype data from reference, target, and source populations.
        ref_ind_file str: Name of the file containing sample information from reference populations.
        tgt_ind_file str: Name of the file containing sample information from target populations.
        src_ind_file str: Name of the file containing sample information from source populations.
        anc_allele_file str: Name of the file containing ancestral allele information.

    Returns:
        ref_data dict: Genotype data from reference populations.
        ref_samples list: Sample information from reference populations.
        tgt_data dict: Genotype data from target populations.
        tgt_samples list: Sample information from target populations.
        src_data list: Genotype data from source populations.
        src_samples list: Sample information from source populations.
    """
   
    ref_data = ref_samples = tgt_data = tgt_samples = src_data = src_samples = None
    if ref_ind_file != None: 
        ref_samples = parse_ind_file(ref_ind_file)
        ref_data = read_geno_data(vcf, ref_samples, anc_allele_file, True)
        
    if tgt_ind_file != None: 
        tgt_samples = parse_ind_file(tgt_ind_file)
        tgt_data = read_geno_data(vcf, tgt_samples, anc_allele_file, True)

    if src_ind_file != None:
        src_samples = parse_ind_file(src_ind_file)
        src_data = read_geno_data(vcf, src_samples, anc_allele_file, False)

    if (ref_ind_file != None) and (tgt_ind_file != None):
        chr_names = tgt_data.keys()
        for c in chr_names:
            # Remove variants fixed in both the reference and target individuals
            ref_fixed_variants = np.sum(ref_data[c]['GT'].is_hom_alt(),axis=1) == len(ref_samples)
            tgt_fixed_variants = np.sum(tgt_data[c]['GT'].is_hom_alt(),axis=1) == len(tgt_samples)
            fixed_index = np.logical_and(ref_fixed_variants, tgt_fixed_variants)
            index = np.logical_not(fixed_index)
            fixed_pos =ref_data[c]['POS'][fixed_index]
            ref_data = filter_data(ref_data, c, index)
            tgt_data = filter_data(tgt_data, c, index)
            if src_ind_file != None:
                index = np.logical_not(np.in1d(src_data[c]['POS'], fixed_pos))
                src_data = filter_data(src_data, c, index)

    return ref_data, ref_samples, tgt_data, tgt_samples, src_data, src_samples

#@profile
def get_ref_alt_allele(ref, alt, pos):
    """
    Description:
        Helper function to index REF and ALT alleles with genomic positions.

    Arguments:
        ref list: REF alleles.
        alt list: ALT alleles.
        pos list: Genomic positions.

    Returns:
        ref_allele dict: REF alleles.
        alt_allele dict: ALT alleles.
    """
    
    ref_allele = dict()
    alt_allele = dict()

    for i in range(len(pos)):
        r = ref[i]
        a = alt[i]
        p = pos[i]
        ref_allele[p] = r
        alt_allele[p] = a
   
    return ref_allele, alt_allele

#@profile
def read_anc_allele(anc_allele_file):
    """
    Description:
        Helper function to read ancestral allele information from files.

    Arguments:
        anc_allele_file str: Name of the BED file containing ancestral allele information.

    Returns:
        anc_allele dict: Ancestral allele information.
    """

    anc_allele = dict()
    with open(anc_allele_file, 'r') as f:
        for line in f.readlines():
            e = line.rstrip().split()
            if e[0] not in anc_allele: anc_allele[e[0]] = dict()
            anc_allele[e[0]][int(e[2])] = e[3]

    if not anc_allele: raise Exception(f'No ancestral allele is found! Please check your data.')
    
    return anc_allele

#@profile
def check_anc_allele(data, anc_allele, c):
    """
    Description:
        Helper function to check whether the REF or ALT allele is the ancestral allele.
        If the ALT allele is the ancestral allele, then the genotypes in this position will be flipped.
        If neither the REF nor ALT allele is the ancestral allele, then this position will be removed.
        If a position has no the ancestral allele information, the this position will be removed.

    Arguments:
        data dict: Genotype data for checking ancestral allele information.
        anc_allele dict: Ancestral allele information for checking.
        c str: Name of the chromosome.

    Returns:
        data dict: Genotype data after checking.
    """

    ref_allele, alt_allele = get_ref_alt_allele(data[c]['REF'], data[c]['ALT'], data[c]['POS'])
    # Remove variants not in the ancestral allele file
    intersect_snps = np.intersect1d(list(ref_allele.keys()), list(anc_allele[c].keys()))
    # Remove variants that neither the ref allele nor the alt allele is the ancestral allele
    removed_snps = []
    # Flip variants that the alt allele is the ancestral allele
    flipped_snps = []

    for v in intersect_snps:
        if (anc_allele[c][v] != ref_allele[v]) and (anc_allele[c][v] != alt_allele[v]): removed_snps.append(v)
        elif (anc_allele[c][v] == alt_allele[v]): flipped_snps.append(v)

    intersect_snps = np.in1d(data[c]['POS'], intersect_snps)
    data = filter_data(data, c, intersect_snps)

    if len(removed_snps) != 0:
        remained_snps = np.logical_not(np.in1d(data[c]['POS'], removed_snps))
        data = filter_data(data, c, remained_snps)

    is_flipped_snps = np.in1d(data[c]['POS'], flipped_snps)
    # Assume no missing data
    for i in range(len(data[c]['POS'])):
        if is_flipped_snps[i]:
            data[c]['GT'][i] = allel.GenotypeVector(abs(data[c]['GT'][i]-1))

    return data

#@profile
def read_mapped_region_file(mapped_region):
    """
    Description:
        Helper function for reading mapped regions from a BED file.

    Arguments:
        mapped_region str: BED file containing mapped regions.

    Returns:
        mapped_intervals dict: Dictionary of tuples containing mapped regions.
    """

    if mapped_region != None:
        mapped_intervals = dict()
        with open(mapped_region, 'r') as m:
            for line in m.readlines():
                elements = line.rstrip().split("\t")
                if elements[0] not in mapped_intervals.keys(): mapped_intervals[elements[0]] = []
                mapped_intervals[elements[0]].append((int(elements[1]), int(elements[2])))
    else: mapped_intervals = None

    return mapped_intervals

#@profile
def cal_matchpct(chr_name, mapped_intervals, data, src_data, tgt_ind_index, src_ind_index, hap_index, win_start, win_end, sample_size):
    """
    Description:
        Helper function to calculate match percents in a given individual.

    Arguments:
        chr_name str: Name of the chromosomes.
        mapped_intervals: Dictionary of tuples containing mapped regions across the genome.
        data dict: Genotype data from individuals to be calculated match percents.
        src_data dict: Genotype data from source populations.
        tgt_ind_index int: Index of the target individual for calculating match percent.
        src_ind_index int: Index of the source individual for calculating match percent.
        hap_index int: Index of the haplotype for calculating match percent.
        win_start int: Start position of the window for calculating match percent.
        win_end int: End position of the window for calculating match percent.
        sample_size int: Number of individuals analyzed.

    Returns:
        res list: List containing statistics for the given haplotype.
    """

    if (win_start == 'NA') and (win_end == 'NA'): 
        return ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']

    win_start = int(win_start)
    win_end = int(win_end)

    hap_match_pct = 'NA'
    gt = data[chr_name]['GT']
    pos = data[chr_name]['POS']
    src_gt = src_data[chr_name]['GT']
    src_pos = src_data[chr_name]['POS']
    res = []
    
    sub_snps = np.where((pos>=win_start) & (pos<=win_end))[0]
    sub_gt = gt[sub_snps][:,tgt_ind_index]
    sub_pos = pos[sub_snps]

    sub_src_snps = np.where((src_pos>=win_start) & (src_pos<=win_end))[0]
    sub_src_gt = src_gt[sub_src_snps][:,src_ind_index]
    sub_src_pos = src_pos[sub_src_snps]
    missing_index = sub_src_gt.is_missing()
    sub_src_gt = sub_src_gt[~missing_index]
    sub_src_pos = sub_src_pos[~missing_index]
    src_hom_variants = sub_src_pos[sub_src_gt.is_hom_alt()]
    src_het_variants = sub_src_pos[sub_src_gt.is_het()]
    src_variants = np.concatenate((src_hom_variants, src_het_variants))

    hap = sub_gt[:,hap_index]
    hap_len = win_end - win_start
    hap_mapped_len = _cal_mapped_len(mapped_intervals, chr_name, win_start, win_end)
    hap_variants_num, hap_site_num, hap_match_src_allele_num, hap_sfs, hap_match_pct = _cal_hap_stats(gt[sub_snps], hap, sub_pos, src_variants, src_hom_variants, src_het_variants, sample_size)

    res.append(hap_variants_num)
    res.append(hap_len)
    res.append(hap_mapped_len)
    res.append(hap_match_src_allele_num)
    res.append(hap_site_num)
    res.append(hap_sfs)
    res.append(hap_match_pct)

    return res

#@profile
def _cal_mapped_len(mapped_intervals, chr_name, win_start, win_end):
    """
    Description:
        Helper function for calculating length of mapped region in a given window.

    Arguments:
        mapped_intervals dict: Dictionary of tuples containing mapped regions across the genome.
        chr_name str: Name of the chromosome for calculating the length of the mapped region.
        win_start int: Start position of a window.
        win_end int: End position of a window.

    Returns:
        mapped_len int: Length of the mapped region.
    """

    if (mapped_intervals == None) or (chr_name not in mapped_intervals.keys()):
        mapped_len = win_end - win_start
    else:
        overlaps = [(idx[0], idx[1]) for idx in mapped_intervals[chr_name] if (win_start>=idx[0] and win_start<=idx[1]) or (win_end>=idx[0] and win_end<=idx[1]) or (win_start<=idx[0] and win_end>=idx[1])]
        mapped_len = 0
        for idx in overlaps:
            if (idx[0]<=win_start) and (idx[1]>=win_end): mapped_len += win_end - win_start
            elif (idx[0]>=win_start) and (idx[1]<=win_end): mapped_len += idx[1] - idx[0]
            elif (idx[0]<=win_start) and (idx[1]<=win_end): mapped_len += idx[1] - win_start
            else: mapped_len += win_end - idx[0]

    mapped_len = mapped_len // 1000 * 1000
        
    return mapped_len
     
#@profile
def _cal_hap_stats(gt, hap, pos, src_variants, src_hom_variants, src_het_variants, sample_size):
    """
    Description:
        Helper function for calculating statistics for a haplotype.

    Arguments:
        gt allel.GenotypeArray: Genotype data for all the haplotypes within the same window of the haplotype to be analyzed.
        hap allel.GenotypeVector: Genotype data for the haplotype to be analyzed.
        pos list: List containing positions of variants on the haplotype.
        src_variants list: List containing positions of variants on the individual from the source population.
        src_hom_variants list: List containing positions of homozygous variants on the individual from the source population.
        src_het_variants list: List containing positions of heterozygous variants on the individual from the source population.
        sample_size int: Number of individuals analyzed.

    Returns:
        hap_variants_num int: Number of SNPs with derived alleles on the haplotype.
        hap_site_num int: Number of SNPs with derived alleles either on the haplotype or the source genomes.
        hap_match_src_allele_num int: Number of SNPs with derived alleles both on the haplotype and the source genomes.
        hap_sfs int: Average number of derived variants per site per haplotype.
        hap_match_pct float: Match percent of the haplotype.
        sample_size int: Number of individuals analyzed.
    """

    if hap is None: return 'NA', 'NA', 'NA', 'NA', 'NA'
    else:
        hap_variants = pos[np.equal(hap, 1)]
        hap_variants_num = len(hap_variants)
        # Assume the alternative allele is the derived allele
        hap_shared_src_hom_site_num = len(np.intersect1d(hap_variants, src_hom_variants))
        hap_shared_src_het_site_num = len(np.intersect1d(hap_variants, src_het_variants))
        hap_site_num = len(np.union1d(hap_variants, src_variants))
        hap_match_src_allele_num = hap_shared_src_hom_site_num + 0.5*hap_shared_src_het_site_num
        hap_shared_src_site_num = hap_shared_src_hom_site_num + hap_shared_src_het_site_num
        if hap_site_num != 0: hap_match_pct = round(hap_match_src_allele_num/hap_site_num, 6)
        else: hap_match_pct = 'NA'

        hap_sfs = np.sum(np.sum(gt[hap == 1], axis=2), axis=1)
        if hap_sfs.size != 0:
            hap_sfs_mean = np.mean(hap_sfs)
            # See https://stackoverflow.com/questions/10825926/python-3-x-rounding-behavior
            #if not np.isnan(sfs_mean): sfs_mean = int(round(sfs_mean))
            #if not np.isnan(hap_sfs_mean): hap_sfs = int(int(py2round(hap_sfs_mean))/10*108)
            #if not np.isnan(hap_sfs_mean): hap_sfs = int(py2round(hap_sfs_mean))/(2*sample_size)
            if not np.isnan(hap_sfs_mean): hap_sfs = round(hap_sfs_mean/(2*sample_size), 6)
        else:
            hap_sfs = np.nan

    return hap_variants_num, hap_site_num, hap_match_src_allele_num, hap_sfs, hap_match_pct

def py2round(x, d=0):
    p = 10 ** d
    if x > 0: return float(math.floor((x * p) + 0.5))/p
    elif x < 0: float(math.ceil((x * p) - 0.5))/p
    else: return 0.0
