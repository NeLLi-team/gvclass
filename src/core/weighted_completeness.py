"""
Sophisticated weighted completeness metrics with ML enhancements.

This module implements conservation-weighted completeness scoring using
marker frequency statistics and adaptive weighting strategies.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
from functools import lru_cache
from sklearn.preprocessing import RobustScaler, MinMaxScaler
from sklearn.ensemble import IsolationForest

logger = logging.getLogger(__name__)


class WeightedCompletenessCalculator:
    """
    Sophisticated weighted completeness calculator with ML enhancements.

    Features:
    - Conservation-based weighting using marker frequency statistics
    - Adaptive weighting strategies (linear, exponential, information-theoretic)
    - Smart caching with LRU for performance
    - Outlier detection for quality control
    - Graceful fallbacks when marker stats unavailable
    - Vectorized operations for efficiency
    """

    def __init__(
        self,
        marker_stats_path: Path,
        weight_strategy: str = "information_theoretic",
        enable_adaptive_scaling: bool = True,
        cache_size: int = 128,
    ):
        """
        Initialize weighted completeness calculator.

        Args:
            marker_stats_path: Path to marker_stats.tsv file
            weight_strategy: One of 'linear', 'exponential', 'information_theoretic'
            enable_adaptive_scaling: Whether to use ML-based adaptive weight scaling
            cache_size: LRU cache size for marker statistics
        """
        self.marker_stats_path = marker_stats_path
        self.weight_strategy = weight_strategy
        self.enable_adaptive_scaling = enable_adaptive_scaling
        self.cache_size = cache_size
        self._marker_stats_df = None
        self._scaler = None
        self._outlier_detector = None

        # Initialize ML components if adaptive scaling is enabled
        if self.enable_adaptive_scaling:
            self._initialize_ml_components()

    def _initialize_ml_components(self):
        """Initialize ML components for adaptive scaling."""
        # Use RobustScaler for better handling of outliers in conservation percentages
        self._scaler = RobustScaler()

        # Isolation Forest for detecting unusual marker patterns
        self._outlier_detector = IsolationForest(
            contamination=0.1, random_state=42, n_estimators=100
        )

    @property
    @lru_cache(maxsize=1)
    def marker_stats_df(self) -> Optional[pd.DataFrame]:
        """
        Load and cache marker statistics DataFrame.

        Returns:
            DataFrame with marker statistics or None if file not found
        """
        try:
            if not self.marker_stats_path.exists():
                logger.warning(f"Marker stats file not found: {self.marker_stats_path}")
                return None

            df = pd.read_csv(self.marker_stats_path, sep="\t")

            # Validate required columns
            required_cols = ["marker_name", "order_name", "percent_genomes_with_marker"]
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                logger.error(f"Missing columns in marker stats: {missing_cols}")
                return None

            # Clean and validate data
            df = df.dropna(subset=required_cols)
            df = df[df["percent_genomes_with_marker"] > 0]  # Remove invalid percentages

            logger.info(
                f"Loaded marker stats for {len(df)} markers across {df['order_name'].nunique()} orders"
            )
            logger.info(f"Available orders: {sorted(df['order_name'].unique())}")
            logger.info(f"Sample markers: {df['marker_name'].head(10).tolist()}")
            return df

        except Exception as e:
            logger.error(f"Error loading marker stats: {e}")
            return None

    def _calculate_base_weights(
        self, conservation_percentages: np.ndarray
    ) -> np.ndarray:
        """
        Calculate base weights using the specified strategy.

        Args:
            conservation_percentages: Array of conservation percentages (0-100)

        Returns:
            Array of weights
        """
        # Convert to 0-1 scale
        conservation_rates = conservation_percentages / 100.0

        if self.weight_strategy == "linear":
            weights = conservation_rates

        elif self.weight_strategy == "exponential":
            # Exponential weighting emphasizes highly conserved markers
            # Use temperature parameter σ = 20 for reasonable scaling
            weights = np.exp(conservation_percentages / 20.0)

        elif self.weight_strategy == "information_theoretic":
            # Information-theoretic weighting based on Shannon entropy
            # More conserved markers carry more information
            weights = -np.log(
                1 - conservation_rates + 1e-10
            )  # Add small epsilon to avoid log(0)

        else:
            raise ValueError(f"Unknown weight strategy: {self.weight_strategy}")

        return weights

    def _apply_adaptive_scaling(
        self, weights: np.ndarray, taxonomic_order: str, marker_names: List[str]
    ) -> np.ndarray:
        """
        Apply ML-based adaptive scaling to weights.

        Args:
            weights: Base weights
            taxonomic_order: Taxonomic order for context-aware scaling
            marker_names: List of marker names

        Returns:
            Adaptively scaled weights
        """
        if not self.enable_adaptive_scaling or self.marker_stats_df is None:
            return weights

        try:
            # Get order-specific marker statistics for context
            order_stats = self.marker_stats_df[
                self.marker_stats_df["order_name"] == taxonomic_order
            ]

            if len(order_stats) < 5:  # Need minimum samples for meaningful scaling
                return weights

            # Create feature matrix for adaptive scaling
            order_conservations = order_stats[
                "percent_genomes_with_marker"
            ].values.reshape(-1, 1)

            # Fit scaler on order-specific data
            order_scaler = RobustScaler()
            order_scaler.fit(order_conservations)

            # Scale weights based on order context
            weights_reshaped = weights.reshape(-1, 1)
            scaled_weights = order_scaler.transform(weights_reshaped).flatten()

            # Apply min-max normalization to ensure positive weights
            min_max_scaler = MinMaxScaler(feature_range=(0.1, 2.0))
            final_weights = min_max_scaler.fit_transform(
                scaled_weights.reshape(-1, 1)
            ).flatten()

            logger.debug(
                f"Applied adaptive scaling for {taxonomic_order}: "
                f"weight range {final_weights.min():.3f}-{final_weights.max():.3f}"
            )

            return final_weights

        except Exception as e:
            logger.warning(f"Adaptive scaling failed, using base weights: {e}")
            return weights

    def _detect_outliers(self, marker_data: pd.DataFrame) -> np.ndarray:
        """
        Detect outlier markers using Isolation Forest.

        Args:
            marker_data: DataFrame with marker statistics

        Returns:
            Boolean array indicating outliers (True = outlier)
        """
        if not self.enable_adaptive_scaling or len(marker_data) < 10:
            return np.zeros(len(marker_data), dtype=bool)

        try:
            # Create feature matrix for outlier detection
            features = marker_data[["percent_genomes_with_marker"]].values

            # Add duplication factor if available
            if "avg_duplication_factor" in marker_data.columns:
                dup_features = marker_data[["avg_duplication_factor"]].values
                features = np.hstack([features, dup_features])

            # Detect outliers
            outliers = self._outlier_detector.fit_predict(features) == -1

            if outliers.any():
                logger.debug(f"Detected {outliers.sum()} outlier markers")

            return outliers

        except Exception as e:
            logger.warning(f"Outlier detection failed: {e}")
            return np.zeros(len(marker_data), dtype=bool)

    @lru_cache(maxsize=128)
    def get_marker_weights(
        self, taxonomic_order: str, marker_names_tuple: Tuple[str, ...]
    ) -> Dict[str, float]:
        """
        Get cached marker weights for a specific taxonomic order.

        Args:
            taxonomic_order: Taxonomic order name
            marker_names_tuple: Tuple of marker names (for caching)

        Returns:
            Dictionary mapping marker names to weights
        """
        marker_names = list(marker_names_tuple)

        if self.marker_stats_df is None:
            # Fallback to uniform weights
            logger.debug("Using uniform weights (marker stats unavailable)")
            return {marker: 1.0 for marker in marker_names}

        try:
            # Filter marker stats for this order and requested markers
            order_markers = self.marker_stats_df[
                (self.marker_stats_df["order_name"] == taxonomic_order)
                & (self.marker_stats_df["marker_name"].isin(marker_names))
            ].copy()

            if order_markers.empty:
                logger.info(
                    f"No order-specific stats for '{taxonomic_order}', trying fallback"
                )
                # Try fallback to any order if order-specific data unavailable
                order_markers = self.marker_stats_df[
                    self.marker_stats_df["marker_name"].isin(marker_names)
                ].copy()

                if order_markers.empty:
                    logger.info(
                        f"No marker stats found for {taxonomic_order}, using uniform weights. Requested markers: {marker_names[:5]}..."
                    )
                    return {marker: 1.0 for marker in marker_names}

            # Detect and handle outliers
            outliers = self._detect_outliers(order_markers)
            if outliers.any():
                # Reduce influence of outlier markers
                order_markers.loc[outliers, "percent_genomes_with_marker"] = (
                    order_markers["percent_genomes_with_marker"].median()
                )

            # Calculate base weights
            conservation_percentages = order_markers[
                "percent_genomes_with_marker"
            ].values
            base_weights = self._calculate_base_weights(conservation_percentages)

            # Apply adaptive scaling
            adaptive_weights = self._apply_adaptive_scaling(
                base_weights, taxonomic_order, marker_names
            )

            # Create weight dictionary
            weight_dict = dict(zip(order_markers["marker_name"], adaptive_weights))

            # Add uniform weights for markers not in stats
            for marker in marker_names:
                if marker not in weight_dict:
                    weight_dict[marker] = 1.0

            logger.debug(
                f"Generated weights for {len(weight_dict)} markers in {taxonomic_order}"
            )
            return weight_dict

        except Exception as e:
            logger.error(f"Error calculating marker weights: {e}")
            return {marker: 1.0 for marker in marker_names}

    def calculate_weighted_completeness(
        self,
        marker_counts: Dict[str, int],
        taxonomic_order: str,
        expected_markers: List[str],
    ) -> Tuple[float, float, Dict[str, any]]:
        """
        Calculate weighted completeness with comprehensive metrics.

        Args:
            marker_counts: Dictionary of marker name to count
            taxonomic_order: Taxonomic order for context-aware weighting
            expected_markers: List of expected markers for this order

        Returns:
            Tuple of (weighted_completeness, confidence_score, detailed_metrics)
        """
        if not expected_markers:
            return 0.0, 0.0, {}

        # Get marker weights (cached)
        marker_weights = self.get_marker_weights(
            taxonomic_order, tuple(expected_markers)
        )

        # Vectorized calculation for efficiency
        weights_array = np.array(
            [marker_weights[marker] for marker in expected_markers]
        )
        presence_array = np.array(
            [
                1.0 if marker_counts.get(marker, 0) > 0 else 0.0
                for marker in expected_markers
            ]
        )

        # Calculate weighted completeness
        total_weight = np.sum(weights_array)
        weighted_score = (
            np.sum(weights_array * presence_array) / total_weight
            if total_weight > 0
            else 0.0
        )
        weighted_completeness = weighted_score * 100.0

        # Calculate confidence score based on marker coverage and weight distribution
        n_present = np.sum(presence_array)
        n_total = len(expected_markers)

        # Confidence factors
        coverage_factor = n_present / n_total if n_total > 0 else 0.0
        weight_diversity = (
            1.0 - np.std(weights_array) / np.mean(weights_array)
            if np.mean(weights_array) > 0
            else 0.0
        )
        weight_diversity = max(0.0, min(1.0, weight_diversity))  # Clamp to [0,1]

        confidence_score = (coverage_factor * 0.7 + weight_diversity * 0.3) * 100.0

        # Detailed metrics for analysis
        detailed_metrics = {
            "n_expected_markers": n_total,
            "n_present_markers": int(n_present),
            "total_weight": float(total_weight),
            "average_weight": float(np.mean(weights_array)),
            "weight_std": float(np.std(weights_array)),
            "coverage_factor": float(coverage_factor),
            "weight_diversity": float(weight_diversity),
            "weight_strategy": self.weight_strategy,
            "adaptive_scaling_enabled": self.enable_adaptive_scaling,
        }

        logger.debug(
            f"Weighted completeness for {taxonomic_order}: "
            f"{weighted_completeness:.2f}% (confidence: {confidence_score:.2f}%)"
        )

        return weighted_completeness, confidence_score, detailed_metrics


def create_weighted_calculator(database_path: Path) -> WeightedCompletenessCalculator:
    """
    Factory function to create a weighted completeness calculator.

    Args:
        database_path: Path to database directory containing marker_stats.tsv

    Returns:
        Configured WeightedCompletenessCalculator instance
    """
    marker_stats_path = database_path / "marker_stats.tsv"

    return WeightedCompletenessCalculator(
        marker_stats_path=marker_stats_path,
        weight_strategy="information_theoretic",  # Most sophisticated strategy
        enable_adaptive_scaling=True,
        cache_size=128,
    )
