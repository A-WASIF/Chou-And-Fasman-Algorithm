# Implementation of Chou and Fasman method for secondary structure prediction

import numpy as np
import pandas as pd

def print_sequence_with_guide(sequence, helical_seq):
    for i in range(0, len(sequence), 50):  # Iterate over the sequence in chunks of 50 characters
        chunk = sequence[i:i+50]
        guide_chunk = helical_seq[i:i+50]

        print(chunk)
        for j in range(len(guide_chunk)):
            if guide_chunk[j] == "H" or guide_chunk[j] == "S" or guide_chunk[j] == "C":
                print(guide_chunk[j], end = "")
            else:
                print(" ", end = "")
            
        print("\n")


def chou_fasman(sequence):


    # Define the amino acid propensity values for alpha helix
    amino_acid_propensity_alpha = {
                                'E': 1.53,
                                'A': 1.45,
                                'L': 1.34,
                                'H': 1.24,
                                'M': 1.20,
                                'Q': 1.17,
                                'V': 1.14,
                                'W': 1.14,
                                'F': 1.12,
                                'K': 1.07,
                                'I': 1.00,
                                'D': 0.98,
                                'T': 0.82,
                                'R': 0.79,
                                'S': 0.79,
                                'C': 0.77,
                                'N': 0.73,
                                'Y': 0.61,
                                'P': 0.59,
                                'G': 0.53
                                }
    
    # Define the amino acid propensity values for beta sheet
    amino_acid_propensity_beta = {
                                'M': 1.67,
                                'V': 1.65,
                                'I': 1.60,
                                'C': 1.30,
                                'Y': 1.29,
                                'F': 1.28,
                                'Q': 1.23,
                                'L': 1.22,
                                'T': 1.20,
                                'W': 1.19,
                                'A': 0.97,
                                'R': 0.90,
                                'G': 0.81,
                                'D': 0.80,
                                'K': 0.74,
                                'S': 0.72,
                                'H': 0.71,
                                'N': 0.65,
                                'P': 0.62,
                                'E': 0.26
                                }
    
    # Alpha helix prediction
    possible_helix_site = {}
    for i in range(len(sequence)):
        if i + 5 < len(sequence):
            possible_helix_site[sequence[i : i + 6]] = [i, i + 6]

    # filter the possible sites based on the quantity of more makers of alpha helix
    filtered_alpha_sites = {}

    for seq, (start, end) in possible_helix_site.items():
        alpha_count = 0

        for i in range(len(seq)):
            if i < 6:
                if(amino_acid_propensity_alpha[seq[i]] >= 1.0):
                    alpha_count += 1
                
        if(alpha_count >= 4):
            filtered_alpha_sites[seq] = (start, end)


    # Extend the alpha helix sites
    for seq in filtered_alpha_sites.items():

        # Left extension
        start = filtered_alpha_sites[seq[0]][0]

        while start > 0:
            prev = sequence[start-1]
            new_seq = prev + sequence[start : start + 3]

            # Fine score of the new sequence
            score = 0
            for i in range(len(new_seq)):
                score += amino_acid_propensity_alpha[new_seq[i]]
            
            if score < 4.0:
                filtered_alpha_sites[seq[0]] = (start, filtered_alpha_sites[seq[0]][1])
                break
            else:
                start -= 1

        # Right extension
        end = filtered_alpha_sites[seq[0]][1]

        while end < len(sequence):
            next = sequence[end]
            new_seq = sequence[end - 3 : end] + next

            # Find score of the new sequence
            score = 0
            for i in range(len(new_seq)):
                score += amino_acid_propensity_alpha[new_seq[i]]
            
            if score < 4.0:
                filtered_alpha_sites[seq[0]] = (filtered_alpha_sites[seq[0]][0], end)
                break
            else:
                end += 1


    # Assigning the helical regions
    helical_regions = [0] * len(sequence)
    
    for seq in filtered_alpha_sites.keys():
        start = filtered_alpha_sites[seq][0]
        end = filtered_alpha_sites[seq][1]
        for i in range(start, end):
            helical_regions[i] = "H"

    # print Alpha Helix regions
    print("------------------------------------------------- Helical Regions ------------------------------------------------------------")
    print_sequence_with_guide(sequence, helical_regions)
    print("\n\n\n")



    # Beta sheet prediction

    possible_sheet_site = {}
    for i in range(len(sequence)):
        if i + 4 < len(sequence):
            possible_sheet_site[sequence[i : i + 5]] = [i, i + 5]
        
    # filter the possible sites based on the quantity of more makers of beta sheet
    filtered_sheet_sites = {}

    for seq, (start, end) in possible_sheet_site.items():
        sheet_count = 0

        for i in range(len(seq)):
            if i < 5:
                if(amino_acid_propensity_beta[seq[i]] >= 1.0):
                    sheet_count += 1
                
        if(sheet_count >= 3):
            filtered_sheet_sites[seq] = (start, end)
        

    # Extend the beta sheet sites
    for seq in filtered_sheet_sites.items():

        # Left extension
        start = filtered_sheet_sites[seq[0]][0]

        while start > 0:
            prev = sequence[start-1]
            new_seq = prev + sequence[start : start + 3]

            # Find score of the new sequence
            score = 0
            for i in range(len(new_seq)):
                score += amino_acid_propensity_beta[new_seq[i]]
            
            if score < 4.0:
                filtered_sheet_sites[seq[0]] = (start, filtered_sheet_sites[seq[0]][1])
                break
            else:
                start -= 1


        # Right extension
        end = filtered_sheet_sites[seq[0]][1]

        while end < len(sequence):
            next = sequence[end]
            new_seq = sequence[end - 3 : end] + next

            # Fine score of the new sequence
            score = 0
            for i in range(len(new_seq)):
                score += amino_acid_propensity_beta[new_seq[i]]
            
            if score < 4.0:
                filtered_sheet_sites[seq[0]] = (filtered_sheet_sites[seq[0]][0], end)
                break
            else:
                end += 1


    # Assigning the beta sheet regions
    sheet_regions = [0] * len(sequence)

    for seq in filtered_sheet_sites.keys():
        start = filtered_sheet_sites[seq][0]
        end = filtered_sheet_sites[seq][1]
        for i in range(start, end):
            sheet_regions[i] = "S"

    # print beta sheet regions
    print("------------------------------------------------- Beta Sheet Regions ------------------------------------------------------------")
    print_sequence_with_guide(sequence, sheet_regions)
    print("\n\n\n")


    # Conflicting regions
    conflicting_regions = [0] * len(sequence)

    for i in range(len(sequence)):
        if helical_regions[i] == "H" and sheet_regions[i] == "S":
            conflicting_regions[i] = "C"

    # print conflicting regions
    print("------------------------------------------------- Conflicting Regions ------------------------------------------------------------")
    print_sequence_with_guide(sequence, conflicting_regions)
    print("\n\n\n")


    # Conflict resolution between Alpha helix and Beta sheet
    resolved_regions = [0] * len(sequence)

    i = 0
    while i < len(sequence):
        if conflicting_regions[i] == "C":
            j = i
            while j < len(sequence) and conflicting_regions[j] == "C":
                j += 1
            
            helix_score = 0
            sheet_score = 0

            for k in range(i, j):
                helix_score += amino_acid_propensity_alpha[sequence[k]]
                sheet_score += amino_acid_propensity_beta[sequence[k]]
            
            if helix_score > sheet_score:
                for k in range(i, j):
                    resolved_regions[k] = "H"
            
            else:
                for k in range(i, j):
                    resolved_regions[k] = "S"

            i = j
        
        else:
            if(helical_regions[i] == "H"):
                resolved_regions[i] = "H"
            elif(sheet_regions[i] == "S"):
                resolved_regions[i] = "S"

            i += 1


    print("------------------------------------------------- Final Resolved Regions ------------------------------------------------------------")
    print_sequence_with_guide(sequence, resolved_regions)
    print("\n\n\n")




if __name__ == "__main__":
    sequence = str("MNASSEGESFAGSVQIPGGTTVLVELTPDIHICGICKQQFNNLDAFVAHKQSGCQLTG" +
                   "TSAAAPSTVQFVSEETVPATQTQTTTRTITSETQTITVSAPEFVFEHGYQTYLPTESN" +
                   "ENQTATVISLPAKSRTKKPTTPPAQKRLNCCYPGCQFKTAYGMKDMERHLKIHTGDKP" +
                   "HKCEVCGKCFSRKDKLKTHMRCHTGVKPYKCKTCDYAAADSSSLNKHLRIHSDERPFK" +
                   "CQICPYASRNSSQLTVHLRSHTASELDDDVPKANCLSTESTDTPKAPVITLPSEAREQ" +
                   "MATLGERTFNCCYPGCHFKTVHGMKDLDRHLRIHTGDKPHKCEFCDKCFSRKDNLTMH" +
                   "MRCHTSVKPHKCHLCDYAAVDSSSLKKHLRIHSDERPYKCQLCPYASRNSSQLTVHLR" +
                   "SHTGDTPFQCWLCSAKFKISSDLKRHMIVHSGEKPFKCEFCDVRCTMKANLKSHIRIK" +
                   "HTFKCLHCAFQGRDRADLLEHSRLHQADHPEKCPECSYSCSSAAALRVHSRVHCKDRP" +
                   "FKCDFCSFDTKRPSSLAKHVDKVHRDEAKTENRAPLGKEGLREGSSQHVAKIVTQRAF" +
                   "RCETCGASFVRDDSLRCHKKQHSDQSENKNSDLVTFPPESGASGQLSTLVSVGQLEAPL" +
                   "EPSQDL")
    
    chou_fasman(sequence)
