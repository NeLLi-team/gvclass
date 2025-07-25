"""
Tree analysis module for extracting nearest neighbors from phylogenetic trees.
"""

import logging
from pathlib import Path
from typing import Dict, Tuple, Optional
from collections import Counter
from ete3 import Tree

logger = logging.getLogger(__name__)


class TreeAnalyzer:
    """Analyze phylogenetic trees to find nearest neighbors."""

    def __init__(self, labels_file: Path):
        """
        Initialize tree analyzer with taxonomy labels.

        Args:
            labels_file: Path to labels file mapping genome IDs to taxonomy
        """
        self.labels_file = labels_file
        self.labels_dict = self.load_labels()

    def load_labels(self) -> Dict[str, str]:
        """Load taxonomy labels from file."""
        labels = {}
        try:
            with open(self.labels_file, "r") as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        genome_id = parts[0]
                        # Store the full taxonomy string
                        labels[genome_id] = parts[1]
        except Exception as e:
            logger.error(f"Error loading labels: {e}")

        return labels

    def get_neighbors_from_tree(
        self, tree_file: Path, query_id: str
    ) -> Dict[str, Dict[str, float]]:
        """
        Extract nearest neighbors for query sequences in a tree.

        Args:
            tree_file: Path to newick tree file
            query_id: Query genome identifier

        Returns:
            Dictionary mapping query proteins to their nearest neighbors and distances
        """
        neighbors = {}

        try:
            tree = Tree(str(tree_file))

            # Find all query nodes in the tree
            query_nodes = []
            for node in tree.traverse():
                if node.is_leaf() and node.name.startswith(query_id):
                    query_nodes.append(node)

            # For each query node, find its nearest non-query neighbor
            for query_node in query_nodes:
                nearest, distance = self._find_nearest_neighbor(
                    tree, query_node, query_id
                )
                if nearest:
                    neighbors[query_node.name] = {nearest: distance}

        except Exception as e:
            logger.error(f"Error processing tree {tree_file}: {e}")

        return neighbors

    def _find_nearest_neighbor(
        self, tree: Tree, query_node: Tree, query_id: str
    ) -> Tuple[Optional[str], float]:
        """
        Find the nearest non-query neighbor for a query node.

        Args:
            tree: The phylogenetic tree
            query_node: The query node
            query_id: Query genome identifier

        Returns:
            Tuple of (neighbor_name, distance) or (None, inf)
        """
        min_distance = float("inf")
        nearest_neighbor = None

        # Get all non-query leaf nodes
        for node in tree.iter_leaves():
            if not node.name.startswith(query_id):
                try:
                    distance = tree.get_distance(query_node.name, node.name)
                    if distance < min_distance:
                        min_distance = distance
                        nearest_neighbor = node.name
                except Exception:
                    # Skip if distance calculation fails
                    continue

        return nearest_neighbor, min_distance

    def get_taxonomy_consensus(
        self,
        all_neighbors: Dict[str, Dict[str, Dict[str, float]]],
        tax_level: str = "order",
    ) -> Optional[str]:
        """
        Get consensus taxonomy from all nearest neighbors across markers.

        Args:
            all_neighbors: Dictionary mapping markers to their neighbors
            tax_level: Taxonomic level ('order', 'family', or 'genus')

        Returns:
            Consensus taxonomy or None
        """

        # Collect all taxonomies
        taxonomies = []

        for marker, query_neighbors in all_neighbors.items():
            for query_protein, neighbors in query_neighbors.items():
                for neighbor, distance in neighbors.items():
                    # Extract genome ID from protein ID
                    genome_id = neighbor.split("|")[0]

                    if genome_id in self.labels_dict:
                        tax_parts = self.labels_dict[genome_id].split("|")
                        # tax_parts: [Kingdom, Class, Order, Family, Genus, Species]
                        # Map to correct index: kingdom=0, order=2, family=3, genus=4
                        tax_idx_map = {
                            "kingdom": 0,
                            "order": 2,
                            "family": 3,
                            "genus": 4,
                        }
                        tax_idx = tax_idx_map.get(tax_level, 2)

                        if len(tax_parts) > tax_idx:
                            taxonomies.append(tax_parts[tax_idx])

        if not taxonomies:
            return None

        # Get most common taxonomy
        tax_counter = Counter(taxonomies)
        most_common = tax_counter.most_common(1)

        if most_common:
            consensus_tax, count = most_common[0]
            # Require at least 50% agreement
            if count >= len(taxonomies) * 0.5:
                return consensus_tax

        return None

    def process_marker_trees(
        self, tree_dir: Path, query_id: str, output_file: Path
    ) -> Dict[str, Dict[str, Dict[str, float]]]:
        """
        Process all marker trees for a query and write results.

        Args:
            tree_dir: Directory containing tree files
            query_id: Query genome identifier
            output_file: Path to output file

        Returns:
            Dictionary of all neighbors by marker
        """
        all_neighbors = {}

        logger.info(f"Processing marker trees in {tree_dir} for query {query_id}")
        logger.info(f"Output file will be: {output_file}")

        # Process each tree file
        tree_files = list(tree_dir.glob("*.treefile"))
        logger.info(f"Found {len(tree_files)} tree files to process")

        for tree_file in tree_files:
            marker = tree_file.stem
            logger.debug(f"Processing tree for marker {marker}")
            neighbors = self.get_neighbors_from_tree(tree_file, query_id)

            if neighbors:
                all_neighbors[marker] = neighbors
                logger.debug(f"Found {len(neighbors)} neighbors for marker {marker}")
            else:
                logger.debug(f"No neighbors found for marker {marker}")

        logger.info(f"Total markers with neighbors: {len(all_neighbors)}")

        # Write results
        if all_neighbors:
            logger.info(f"Writing tree_nn results to {output_file}")
            with open(output_file, "w") as f:
                f.write("# Tree-based nearest neighbors\n")
                f.write("marker\tquery_protein\tnearest_neighbor\tdistance\ttaxonomy\n")

                for marker, query_neighbors in all_neighbors.items():
                    for query_protein, neighbors in query_neighbors.items():
                        for neighbor, distance in neighbors.items():
                            genome_id = neighbor.split("|")[0]
                            tax_str = "unknown"

                            if genome_id in self.labels_dict:
                                tax_parts = self.labels_dict[genome_id].split("|")
                                # Format: kingdom;order;family;genus
                                if len(tax_parts) >= 5:
                                    tax_str = f"{tax_parts[0]};{tax_parts[2]};{tax_parts[3]};{tax_parts[4]}"

                            f.write(
                                f"{marker}\t{query_protein}\t{neighbor}\t{distance:.6f}\t{tax_str}\n"
                            )
            logger.info("Successfully wrote tree_nn file")
        else:
            logger.warning(
                "No neighbors found for any markers - not writing tree_nn file"
            )

        return all_neighbors
