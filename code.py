#In genomics, N90 is a metric used to assess the quality of an assembly or sequencing data. 
#It represents the contig or read length at which at least 90% of the total assembly length is contained in contigs or reads of that length or longer.
#Similar to N50, the N90 value is calculated by sorting the contig or read lengths in descending order and then summing the lengths cumulatively. 
#The process continues until the cumulative length reaches or exceeds 90% of the total assembly length. The N90 value is the length of the contig or read at this point.


def calculate_n90(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    n90_length = total_length * 0.9

    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= n90_length:
            return length

    return None  # N90 not found


# Example usage
contig_lengths = [500, 1000, 750, 200, 300, 900, 1500, 1200]
n90 = calculate_n90(contig_lengths)
print("N90 value:", n90)


#In genomics, NG50 (N50 based on the reference genome) is a metric used to assess the quality of an assembly in comparison to a reference genome. 
#It represents the contig or scaffold length at which at least 50% of the reference genome is covered by the assembly.
#To calculate NG50, the contig or scaffold lengths are first sorted in descending order. Then, the lengths are summed cumulatively until the cumulative length reaches or exceeds 50% of the reference genome length. 
#The NG50 value is the length of the contig or scaffold at this point.

def calculate_ng50(reference_length, lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    ng50_length = reference_length * 0.5

    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= ng50_length:
            return length

    return None  # NG50 not found


# Example usage
reference_length = 1000000
contig_lengths = [5000, 10000, 7500, 20000, 30000, 9000, 15000, 12000]
ng50 = calculate_ng50(reference_length, contig_lengths)
print("NG50 value:", ng50)



#In the field of genomics, D50 (Duplication 50) is a metric used to evaluate the level of duplication or repetitive sequences within an assembly or sequencing data. 
#It represents the contig or read length at which 50% of the total assembly or sequencing length is contained within duplicated or repetitive regions.
#To calculate the D50 value, the contig or read lengths are sorted in descending order. Then, the lengths are summed cumulatively while considering the duplicated or repetitive regions. 
#The calculation continues until the cumulative length reaches or exceeds 50% of the total assembly or sequencing length. The D50 value is the length of the contig or read at this point.

def calculate_d50(lengths, duplicated_regions):
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    d50_length = total_length * 0.5

    cumulative_length = 0
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= d50_length:
            return length

    return None  # D50 not found


# Example usage
contig_lengths = [500, 1000, 750, 200, 300, 900, 1500, 1200]
duplications = [False, True, False, False, True, True, False, False]
d50 = calculate_d50(contig_lengths, duplications)
print("D50 value:", d50)


#U50 metric represents the length of the smallest contig that contains 50% of the cumulative sum of all unique, target-specific contigs.

def calculate_u50(contig_lengths):
    sorted_lengths = sorted(contig_lengths)
    total_length = sum(sorted_lengths)
    u50_threshold = total_length * 0.5

    cumulative_length = 0
    u50 = None
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= u50_threshold:
            u50 = length
            break

    return u50


# Example usage
contig_lengths = [500, 1000, 750, 200, 300, 900, 1500, 1200]
u50 = calculate_u50(contig_lengths)
print("U50 value:", u50)

#UL50 is the number of contigs whose length sum produces U50

def calculate_ul50(contig_lengths, u50):
    sorted_lengths = sorted(contig_lengths)
    cumulative_length = 0
    ul50 = 0
    for length in sorted_lengths:
        cumulative_length += length
        ul50 += 1
        if cumulative_length >= u50:
            break
    return ul50


# Example usage
contig_lengths = [500, 1000, 750, 200, 300, 900, 1500, 1200]
u50 = 2500
ul50 = calculate_ul50(contig_lengths, u50)
print("UL50 value:", ul50)


#UG50 is the length of the smallest contig such that 50% of the reference genome is contained in unique, target-specific contigs of size UG50 or larger.

def calculate_ug50(contig_lengths, reference_genome_size):
    sorted_lengths = sorted(contig_lengths)
    total_length = sum(sorted_lengths)
    ug50_threshold = reference_genome_size * 0.5

    cumulative_length = 0
    ug50 = None
    for length in sorted_lengths:
        cumulative_length += length
        if cumulative_length >= ug50_threshold:
            ug50 = length
            break

    return ug50


# Example usage
contig_lengths = [500, 1000, 750, 200, 300, 900, 1500, 1200]
reference_genome_size = 3000
ug50 = calculate_ug50(contig_lengths, reference_genome_size)
print("UG50 value:", ug50)


#UG50% is the estimated percent coverage length of the UG50 in direct relation to the length of the reference genome. 
#The calculation is (100 × (UG50/Length of reference genome). 
#The UG50%, as a percentage-based metric, can be used to compare assembly results from different samples or studies.

def calculate_ug50_percentage(ug50, reference_genome_length):
    ug50_percentage = (ug50 / reference_genome_length) * 100
    return ug50_percentage


# Example usage
ug50 = 2500
reference_genome_length = 5000
ug50_percentage = calculate_ug50_percentage(ug50, reference_genome_length)
print("UG50% value:", ug50_percentage)

