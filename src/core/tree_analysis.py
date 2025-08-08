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
                line_count = 0
                for line in f:
                    line_count += 1
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        genome_id = parts[0]
                        # Store the full taxonomy string
                        labels[genome_id] = parts[1]
            logger.info(
                f"Loaded {len(labels)} labels from {line_count} lines in {self.labels_file}"
            )
            if len(labels) == 0:
                logger.warning(f"No labels were loaded from {self.labels_file}!")
                # Show first few lines for debugging
                with open(self.labels_file, "r") as f:
                    first_lines = [f.readline() for _ in range(5)]
                logger.warning(f"First 5 lines of labels file: {first_lines}")
        except Exception as e:
            logger.error(f"Error loading labels: {e}")
            logger.error(f"Labels file path: {self.labels_file}")
            logger.error(f"File exists: {self.labels_file.exists()}")

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

        # Check if tree file exists and is not a dummy/skipped tree
        if not tree_file.exists():
            logger.debug(
                f"Tree file {tree_file} does not exist - likely skipped due to insufficient sequences"
            )
            return neighbors

        # Check for skip marker file
        skip_file = tree_file.with_suffix(".SKIPPED")
        if skip_file.exists():
            logger.debug(
                f"Tree building was skipped for {tree_file.stem} - insufficient sequences"
            )
            return neighbors

        try:
            tree = Tree(str(tree_file))

            # Debug: Log tree parsing
            leaf_count = len([n for n in tree.iter_leaves()])
            logger.debug(f"Loaded tree from {tree_file.name} with {leaf_count} leaves")

            # Find all query nodes in the tree
            query_nodes = []
            for node in tree.traverse():
                if node.is_leaf() and node.name.startswith(query_id):
                    query_nodes.append(node)

            if len(query_nodes) == 0:
                # Debug: Show what leaf names we have
                sample_leaves = [n.name for n in tree.iter_leaves()][:5]
                logger.debug(f"No query nodes found for {query_id} in {tree_file.name}")
                logger.debug(f"Sample leaf names: {sample_leaves}")
                logger.debug(f"Looking for nodes starting with: {query_id}")

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
                    # Extract genome ID from protein ID (handle missing | gracefully)
                    if "|" in neighbor:
                        genome_id = neighbor.split("|")[0]
                    else:
                        genome_id = neighbor

                    if genome_id in self.labels_dict:
                        tax_str = self.labels_dict[genome_id]
                        # Split taxonomy string, handling various separators
                        if "|" in tax_str:
                            tax_parts = tax_str.split("|")
                        else:
                            # Handle single value or unknown format
                            tax_parts = [tax_str]

                        # tax_parts: [Kingdom, Class, Order, Family, Genus, Species]
                        # Map to correct index: kingdom=0, order=2, family=3, genus=4
                        tax_idx_map = {
                            "kingdom": 0,
                            "order": 2,
                            "family": 3,
                            "genus": 4,
                        }
                        tax_idx = tax_idx_map.get(tax_level, 2)

                        if len(tax_parts) > tax_idx and tax_parts[tax_idx].strip():
                            taxonomies.append(tax_parts[tax_idx].strip())

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
        skipped_files = list(tree_dir.glob("*.SKIPPED"))
        logger.info(f"Found {len(tree_files)} tree files to process")
        if skipped_files:
            logger.info(
                f"Note: {len(skipped_files)} markers had tree building skipped due to insufficient sequences"
            )

        trees_processed = 0
        trees_skipped = 0

        for tree_file in tree_files:
            marker = tree_file.stem
            logger.debug(f"Processing tree for marker {marker}")

            # Check if this tree was skipped
            if tree_file.with_suffix(".SKIPPED").exists():
                trees_skipped += 1
                logger.debug(
                    f"Skipping {marker} - insufficient sequences for tree building"
                )
                continue

            neighbors = self.get_neighbors_from_tree(tree_file, query_id)

            if neighbors:
                all_neighbors[marker] = neighbors
                logger.debug(f"Found {len(neighbors)} neighbors for marker {marker}")
                trees_processed += 1
            else:
                logger.debug(f"No neighbors found for marker {marker}")

        logger.info(f"Total markers with neighbors: {len(all_neighbors)}")
        logger.info(
            f"Trees processed: {trees_processed}, Trees skipped: {trees_skipped}"
        )

        # Write results
        if all_neighbors:
            logger.info(f"Writing tree_nn results to {output_file}")
            with open(output_file, "w") as f:
                f.write("# Tree-based nearest neighbors\n")
                f.write("marker\tquery_protein\tnearest_neighbor\tdistance\ttaxonomy\n")

                for marker, query_neighbors in all_neighbors.items():
                    for query_protein, neighbors in query_neighbors.items():
                        for neighbor, distance in neighbors.items():
                            # Extract genome ID (handle missing | gracefully)
                            if "|" in neighbor:
                                genome_id = neighbor.split("|")[0]
                            else:
                                genome_id = neighbor
                            tax_str = "unknown"

                            if genome_id in self.labels_dict:
                                label_str = self.labels_dict[genome_id]
                                if "|" in label_str:
                                    tax_parts = label_str.split("|")
                                    # Format: kingdom;order;family;genus (skip empty values)
                                    if len(tax_parts) >= 5:
                                        kingdom = tax_parts[0].strip() or "unknown"
                                        order = (
                                            tax_parts[2].strip()
                                            if len(tax_parts) > 2
                                            else ""
                                        )
                                        family = (
                                            tax_parts[3].strip()
                                            if len(tax_parts) > 3
                                            else ""
                                        )
                                        genus = (
                                            tax_parts[4].strip()
                                            if len(tax_parts) > 4
                                            else ""
                                        )
                                        tax_str = f"{kingdom};{order};{family};{genus}"
                                else:
                                    # Handle non-standard format
                                    tax_str = label_str.strip() or "unknown"

                            f.write(
                                f"{marker}\t{query_protein}\t{neighbor}\t{distance:.6f}\t{tax_str}\n"
                            )
            logger.info("Successfully wrote tree_nn file")
        else:
            logger.warning(
                "No neighbors found for any markers - not writing tree_nn file"
            )

        return all_neighbors
