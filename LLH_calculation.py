from build_tree import tree

def calculate_subtree(subtree_id):
    subtree = tree.subtree(subtree_id)
    return subtree

a = tree.get_node(1023)

print(tree.get_node(1023))

calculate_subtree(1023).show() ##A problem here - all good - 1023 is identifier - not tag




