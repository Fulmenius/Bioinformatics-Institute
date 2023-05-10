# -*- coding: utf-8 -*-
"""
Accessing http://hollywood.mit.edu/GENSCAN.html to obtain exons and introns 
from a nucleotide sequence. 
"""

from bs4 import BeautifulSoup
import requests
import re
import pandas as pd
from typing import IO, Union, Literal
from dataclasses import dataclass

@dataclass
class GenscanOutput:
    status: str
    cds_list: str
    intron_list: pd.DataFrame
    exon_list: pd.DataFrame

def run_genscan(sequence: str = None,
                sequence_file: Union[str, IO] = None,
                organism: Literal["Vertebrate", "Arabidopsis", "Maze"] = "Vertebrate",
                exon_cutoff: float = 1.00,
                sequence_name: str = "") -> GenscanOutput:
                
        FORM_URL = "http://hollywood.mit.edu/cgi-bin/genscanw_py.cgi"
        
        # Choose if the form is to be filled using a file or "manual" input
        if sequence_file:
            data = {"-o":organism,
                    "-e":exon_cutoff,
                    "-n":sequence_name,
                    "-p":"Predicted peptides only",
                    "-u":sequence_file}
        else:
            data = {"-o":organism,
                    "-e":exon_cutoff,
                    "-n":sequence_name,
                    "-p":"Predicted peptides only",
                    "-s":sequence}

        response = requests.post(url=FORM_URL, data=data)
        status = response.status_code
        
        soup = BeautifulSoup(response.content, "lxml")
        text = soup.find_all("pre").__repr__()
        
        # A tiny itsy-bitsy little humble bit of GPT-4 diabolic sorcery below
        
        # Extract exons and introns
        exon_intron_pattern = r"(\d+\.\d+)\s+(Intr|Term|PlyA)\s+(\+|\-)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)?\s+(\S+)?\s+(\S+)?"
        matches = re.findall(exon_intron_pattern, text)

        exon_data = []
        intron_data = []

        for match in matches:
            match = tuple(x if x is not None else '' for x in match)
            if match[1] == "Intr":
                intron_data.append(match)
            elif match[1] in ["Term", "PlyA"]:
                exon_data.append(match[:-1] + ('',))

        # Make sure each row has the correct length
        exon_data = [row if len(row) == 12 else row + ('',) for row in exon_data]
        intron_data = [row if len(row) == 13 else row + ('',) for row in intron_data]

        columns = ["ID", "Type", "Strand", "Begin", "End", "Length", "Frame", "Phase", "I/Ac", "Do/T", "CodRg", "P", "Tscr"]
        exon_df = pd.DataFrame(exon_data, columns=columns[:-1])
        intron_df = pd.DataFrame(intron_data, columns=columns)

        # Extract amino acid sequence
        amino_acid_pattern = r"(?<=predicted_peptide_1\|)(\d+)_aa\n\n(.*\n)*"
        amino_acid_match = re.search(amino_acid_pattern, text, flags=re.MULTILINE)
        amino_acid_sequence = amino_acid_match.group(0).replace("\n", "") if amino_acid_match else ""

        exon_list = exon_df
        intron_list = intron_df
        cds_list = amino_acid_sequence.split("_aa")[1] if len(amino_acid_sequence) > 0 else amino_acid_sequence

        return GenscanOutput(status=status, cds_list=cds_list, intron_list=intron_list, exon_list=exon_list)