from ete3 import Tree

def split_newick_trees(file_path):
    """Split concatenated Newick trees by semicolon."""
    with open(file_path) as f:
        content = f.read().strip()

    trees = []
    current = ""

    for char in content:
        current += char
        if char == ';':
            trees.append(current.strip())
            current = ""

    return trees


def validate_trees(file_path):
    trees = split_newick_trees(file_path)

    print(f"\nTotal trees found: {len(trees)}")

    valid_count = 0
    invalid_count = 0

    reference_taxa = None

    for i, t_str in enumerate(trees):
        try:
            t = Tree(t_str, format=1)

            # --- Check duplicate taxa ---
            leaf_names = t.get_leaf_names()
            if len(leaf_names) != len(set(leaf_names)):
                raise ValueError("Duplicate taxa")

            # --- Check for polytomies ---
            for node in t.traverse():
                if len(node.children) > 2:
                    raise ValueError("Polytomy detected")

            # --- Check taxa consistency ---
            taxa_set = set(leaf_names)
            if reference_taxa is None:
                reference_taxa = taxa_set
            elif taxa_set != reference_taxa:
                raise ValueError("Taxa mismatch")

            valid_count += 1

        except Exception as e:
            print(f"Tree {i+1} INVALID -> {str(e)}")
            invalid_count += 1

    print("\nSUMMARY")
    print(f"Valid trees: {valid_count}")
    print(f"Invalid trees: {invalid_count}")


# ✅ Run it
validate_trees("concat_for_quibl.treefile")