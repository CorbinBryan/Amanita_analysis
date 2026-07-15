#!/usr/bin/env python3

import os

input_file = "/Users/corbinbryan/Desktop/Amanita_analysis/buscos/concat_for_quibl.tree"
output_file = "/Users/corbinbryan/Desktop/Amanita_analysis/buscos/renamed_concat_for_quibl.tree"

def rename_label(label):
    return label.split("_")[0]

def process_newick_line(line):
    result = []
    token = ""

    for char in line:
        if char in "(),:;":
            if token:
                # Only rename tip labels (not branch lengths)
                new_token = rename_label(token)
                result.append(new_token)
                token = ""
            result.append(char)
        else:
            token += char

    if token:
        result.append(rename_label(token))

    return "".join(result)

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        line = line.strip()
        if line:
            new_line = process_newick_line(line)
            outfile.write(new_line + "\n")

print("Done!")
