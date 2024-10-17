import sys
from collections import defaultdict

def read_vcf(vcf_file):
    vcf_dict = defaultdict(list)
    with open(vcf_file, 'r') as vcf:
        next(vcf)
        for line in vcf:
            parts = line.strip().split()
            if parts:
                key = (parts[0], parts[-1])
                vcf_dict[key].append(parts)
                
    return vcf_dict

def read_pileup(pileup_file):
    pileup_dict = {}
    with open(pileup_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 8:  
                continue
            
            chrom = parts[0]
            pos = parts[1]
            # Collect base and read count information
            read_info = [int(parts[3]), parts[4].upper(), int(parts[6]), parts[7].upper()]
            pileup_dict[(chrom, pos)] = read_info
            
    return pileup_dict 

def count_matches(pileup_dict, vcf_dict):
    num = 0
    denom = 0
    for k,v in vcf_dict.items():
        if len(v) < 2:
            continue
        snps = {}
        for snp in v:
            #denom+=1
            chrom_pos = (snp[0], snp[1])
            ref, alt = snp[2], snp[3]
            alleles = [ref,alt]
            GT = snp[5]
            PS = snp[6]
            # If snp is present in pileup (assembly)
            if chrom_pos in pileup_dict:
                read_count1, base1, read_count2, base2 = pileup_dict[chrom_pos]
                bases = [base1, base2]
                
                if read_count1 == 0 or read_count2 == 0:
                    continue
                else:
                    denom+=1
                
                #print(snp, ref, alt, GT, read_count1, base1, read_count2, base2)
                # Check matching conditions
                if read_count1 > 1 and read_count2 > 1:
                    # Both files have multiple reads, check if bases cover the alleles
                    if (ref in base1 and alt in base2) or (ref in base2 and alt in base1):
                        if GT == '0|1':
                            if (ref in base1 and alt in base2) and (ref not in base2 and alt not in base1):
                                assembly_GT = 1
                            elif (ref in base2 and alt in base1) and (ref not in base1 and alt not in base2):
                                assembly_GT = 2
                            else:
                                assembly_GT = 0
                        elif GT == '1|0':
                            if (ref in base1 and alt in base2) and (ref not in base2 and alt not in base1):
                                assembly_GT = 2
                            elif (ref in base2 and alt in base1) and (ref not in base1 and alt not in base2):
                                assembly_GT = 1
                            else:
                                assembly_GT = 0

                        snps[chrom_pos] = assembly_GT
                elif (read_count1 == 1 and read_count2 > 1) or (read_count1 > 1 and read_count2 == 1):
                    if read_count1 == 1 and read_count2 > 1:
                        if GT == '0|1':
                            if (ref == base1 and alt in base2):
                                assembly_GT = 1
                            elif (alt == base1 and ref in base2):
                                assembly_GT = 2
                            else:
                                assembly_GT = 0
                        elif GT == '1|0':
                            if (ref == base1 and alt in base2):
                                assembly_GT = 2
                            elif (alt == base1 and ref in base2):
                                assembly_GT = 1
                            else:
                                assembly_GT = 0

                    elif read_count2 == 1 and read_count1 > 1:
                        if GT == '0|1':
                            if (ref == base2 and alt in base1):
                                assembly_GT = 2
                            elif (alt == base2 and ref in base1):
                                assembly_GT = 1
                            else:
                                assembly_GT = 0
                        elif GT == '1|0':
                            if (ref == base2 and alt in base1):
                                assembly_GT = 1
                            elif (alt == base2 and ref in base1):
                                assembly_GT = 2
                            else:
                                assembly_GT = 0

                    snps[chrom_pos] = assembly_GT   
                elif read_count1 == 1 and read_count2 == 1:
                    if GT == '0|1':
                        if (ref == base1 and alt == base2):
                            assembly_GT = 1
                        elif (alt == base1 and ref == base2):
                            assembly_GT = 2
                        else:
                            assembly_GT = 0
                    elif GT == '1|0':
                        if (ref == base1 and alt == base2):
                            assembly_GT = 2
                        elif (alt == base1 and ref == base2):
                            assembly_GT = 1
                        else:
                            assembly_GT = 0

                    snps[chrom_pos] = assembly_GT
                #print(bases, alleles, GT, assembly_GT)
        # Checking adjacent snps
        keys = list(snps.keys())

        for index in range(len(keys)):
            prev_value = None
            next_value = None
            if index > 0:
                prev_key = keys[index - 1]
                prev_value = snps[prev_key]
            if index < len(keys) - 1:
                next_key = keys[index + 1]
                next_value = snps[next_key]

            current_value = snps[keys[index]]
            if current_value in [prev_value, next_value] and current_value != 0:
                num += 1

    return num, num/denom

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py hetsnp_vcf_file.txt pileup_file.txt ")
        sys.exit(1)

    vcf_file = sys.argv[1]
    pileup_file = sys.argv[2]

    vcf_dict = read_vcf(vcf_file)
    pileup_dict = read_pileup(pileup_file)
    matches, percent = count_matches(pileup_dict, vcf_dict)
    percentage=percent*100
    print(f"Matching bases count: {matches} ({percentage:.2f}%)")

if __name__ == "__main__":
    main()