try:
    from ete3 import Tree
    _HAVE_ETE3 = True
except Exception:
    _HAVE_ETE3 = False
    import re

input_file = "concat_for_quibl.treefile"
output_file = "filtered_trees_for_quibl.treefile"
target_taxa = 8

def _count_leaves_fallback(newick_line):
    # Fallback heuristic: count label occurrences before a ':' (label:branch_length)
    # This works for typical leaf labels like 'Taxon:0.001'.
    return len(re.findall(r"([A-Za-z0-9_.]+):", newick_line))

with open(input_file) as f, open(output_file, "w") as out:
    for line in f:
        line = line.strip()
        if not line:
            continue

        if _HAVE_ETE3:
            t = Tree(line)
            n_leaves = len(t.get_leaves())
        else:
            n_leaves = _count_leaves_fallback(line)

        if n_leaves == target_taxa:
            out.write(line + "\n")