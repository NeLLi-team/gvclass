"""
Comprehensive tests for weighted completeness metrics.

Tests the ML-enhanced weighted completeness calculator with various scenarios
including edge cases, fallback behaviors, and performance characteristics.
"""

import pytest
import numpy as np
import pandas as pd
from pathlib import Path
import tempfile
import logging
from unittest.mock import patch

from src.core.weighted_completeness import (
    WeightedCompletenessCalculator,
    create_weighted_calculator,
)
from src.core.summarize_full import FullSummarizer


@pytest.fixture
def sample_marker_stats():
    """Create sample marker statistics data."""
    return pd.DataFrame(
        {
            "marker_name": [
                "OGv21091",
                "OGv21249",
                "OGv21407",
                "OGv21419",
                "OGv21421",
                "OGv2118",
                "OGv2122",
                "OGv2166",
                "OGv2183",
                "OGv2186",
            ],
            "domain_name": ["MIRUS"] * 5 + ["NCLDV"] * 5,
            "order_name": ["MIRUS"] * 5 + ["Algavirales"] * 5,
            "percent_genomes_with_marker": [
                33.7,
                42.4,
                39.2,
                38.4,
                38.9,
                54.9,
                41.8,
                44.4,
                59.5,
                63.2,
            ],
            "avg_duplication_factor": [
                1.0159,
                1.0116,
                1.0126,
                1.0127,
                1.0063,
                1.0368,
                1.0194,
                1.0380,
                1.0373,
                1.0330,
            ],
        }
    )


@pytest.fixture
def temp_marker_stats_file(sample_marker_stats):
    """Create temporary marker stats file."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        sample_marker_stats.to_csv(f, sep="	", index=False)
        temp_path = Path(f.name)
    yield temp_path
    temp_path.unlink()


@pytest.fixture
def calculator(temp_marker_stats_file):
    """Create WeightedCompletenessCalculator instance."""
    return WeightedCompletenessCalculator(
        marker_stats_path=temp_marker_stats_file,
        weight_strategy="information_theoretic",
        enable_adaptive_scaling=True,
    )

class TestWeightedCompletenessCalculator:
    """Test suite for WeightedCompletenessCalculator class."""

    def test_initialization(self, temp_marker_stats_file):
        """Test calculator initialization with different strategies."""
        strategies = ["linear", "exponential", "information_theoretic"]

        for strategy in strategies:
            calc = WeightedCompletenessCalculator(
                marker_stats_path=temp_marker_stats_file,
                weight_strategy=strategy,
                enable_adaptive_scaling=True,
            )
            assert calc.weight_strategy == strategy
            assert calc.enable_adaptive_scaling is True
            assert calc._scaler is not None
            assert calc._outlier_detector is not None

    def test_marker_stats_loading(self, calculator, sample_marker_stats):
        """Test marker statistics loading and validation."""
        df = calculator.marker_stats_df
        assert df is not None
        assert len(df) == len(sample_marker_stats)
        assert all(
            col in df.columns
            for col in ["marker_name", "order_name", "percent_genomes_with_marker"]
        )

    def test_weight_strategies(self, calculator):
        """Test different weight calculation strategies."""
        conservation_percentages = np.array([25.0, 50.0, 75.0, 90.0])

        # Test linear strategy
        calculator.weight_strategy = "linear"
        linear_weights = calculator._calculate_base_weights(conservation_percentages)
        expected_linear = conservation_percentages / 100.0
        np.testing.assert_array_almost_equal(linear_weights, expected_linear)

        # Test exponential strategy
        calculator.weight_strategy = "exponential"
        exp_weights = calculator._calculate_base_weights(conservation_percentages)
        expected_exp = np.exp(conservation_percentages / 20.0)
        np.testing.assert_array_almost_equal(exp_weights, expected_exp)

        # Test information-theoretic strategy
        calculator.weight_strategy = "information_theoretic"
        info_weights = calculator._calculate_base_weights(conservation_percentages)
        conservation_rates = conservation_percentages / 100.0
        expected_info = -np.log(1 - conservation_rates + 1e-10)
        np.testing.assert_array_almost_equal(info_weights, expected_info)

    def test_weight_caching(self, calculator):
        """Test LRU caching of marker weights."""
        markers = ("OGv21091", "OGv21249", "OGv21407")
        order = "MIRUS"

        # First call should compute weights
        weights1 = calculator.get_marker_weights(order, markers)

        # Second call should use cache
        weights2 = calculator.get_marker_weights(order, markers)

        assert weights1 == weights2
        assert len(weights1) == len(markers)
        assert all(marker in weights1 for marker in markers)

    def test_fallback_behavior(self):
        """Test graceful fallback when marker stats unavailable."""
        # Create calculator with non-existent file
        fake_path = Path("/nonexistent/marker_stats.tsv")
        calculator = WeightedCompletenessCalculator(
            marker_stats_path=fake_path,
            weight_strategy="linear",
            enable_adaptive_scaling=False,
        )

        markers = ("marker1", "marker2", "marker3")
        weights = calculator.get_marker_weights("test_order", markers)

        # Should fallback to uniform weights
        assert all(weight == 1.0 for weight in weights.values())
        assert len(weights) == len(markers)

    def test_outlier_detection(self, calculator, sample_marker_stats):
        """Test outlier detection functionality."""
        # Create data with obvious outliers
        outlier_data = sample_marker_stats.copy()
        outlier_data.loc[0, "percent_genomes_with_marker"] = 99.9  # Outlier
        outlier_data.loc[1, "avg_duplication_factor"] = 10.0  # Outlier

        outliers = calculator._detect_outliers(outlier_data)

        # Should detect some outliers (exact number depends on algorithm)
        assert isinstance(outliers, np.ndarray)
        assert outliers.dtype == bool
        assert len(outliers) == len(outlier_data)

    def test_weighted_completeness_calculation(self, calculator):
        """Test weighted completeness calculation with realistic data."""
        marker_counts = {
            "OGv21091": 1,
            "OGv21249": 1,
            "OGv21407": 0,  # Missing marker
            "OGv21419": 1,
            "OGv21421": 2,  # Duplicated marker
        }

        expected_markers = list(marker_counts.keys())

        weighted_comp, confidence, metrics = calculator.calculate_weighted_completeness(
            marker_counts=marker_counts,
            taxonomic_order="MIRUS",
            expected_markers=expected_markers,
        )

        # Basic validation
        assert 0 <= weighted_comp <= 100
        assert 0 <= confidence <= 100
        assert isinstance(metrics, dict)
        assert "n_expected_markers" in metrics
        assert "n_present_markers" in metrics
        assert metrics["n_expected_markers"] == len(expected_markers)
        assert metrics["n_present_markers"] == 4  # 4 out of 5 markers present

    def test_empty_markers_list(self, calculator):
        """Test behavior with empty markers list."""
        weighted_comp, confidence, metrics = calculator.calculate_weighted_completeness(
            marker_counts={}, taxonomic_order="test", expected_markers=[]
        )

        assert weighted_comp == 0.0
        assert confidence == 0.0
        assert metrics == {}

    def test_adaptive_scaling_disabled(self, temp_marker_stats_file):
        """Test calculator with adaptive scaling disabled."""
        calculator = WeightedCompletenessCalculator(
            marker_stats_path=temp_marker_stats_file,
            weight_strategy="linear",
            enable_adaptive_scaling=False,
        )

        assert calculator._scaler is None
        assert calculator._outlier_detector is None

        # Should still work but without adaptive features
        markers = ("OGv21091", "OGv21249")
        weights = calculator.get_marker_weights("MIRUS", markers)
        assert len(weights) == len(markers)

    def test_unknown_weight_strategy(self, temp_marker_stats_file):
        """Test error handling for unknown weight strategy."""
        calculator = WeightedCompletenessCalculator(
            marker_stats_path=temp_marker_stats_file, weight_strategy="unknown_strategy"
        )

        with pytest.raises(ValueError, match="Unknown weight strategy"):
            calculator._calculate_base_weights(np.array([50.0]))


class TestFullSummarizerIntegration:
    """Test integration of weighted completeness with FullSummarizer."""

    @pytest.fixture
    def temp_database_dir(self):
        """Create temporary database directory with required files."""
        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = Path(temp_dir)

            # Create mock marker stats file
            marker_stats = pd.DataFrame(
                {
                    "marker_name": ["OG21", "OG22", "OG23"],
                    "order_name": ["TestOrder"] * 3,
                    "percent_genomes_with_marker": [50.0, 75.0, 25.0],
                }
            )
            marker_stats.to_csv(db_path / "marker_stats.tsv", sep="\t", index=False)

            # Create mock order completeness file
            order_comp = pd.DataFrame(
                {"Order": ["TestOrder"], "Orthogroups": ["OG21, OG22, OG23"]}
            )
            order_comp.to_csv(db_path / "order_completeness.tab", sep="\t", index=False)

            # Create mock labels file
            labels_content = "genome1\tNCLDV|Class|TestOrder|Family|Genus|Species\n"
            (db_path / "gvclassSeptember25_labels.tsv").write_text(labels_content)

            yield db_path

    def test_summarizer_initialization(self, temp_database_dir):
        """Test FullSummarizer initialization with weighted calculator."""
        summarizer = FullSummarizer(temp_database_dir)

        assert hasattr(summarizer, "weighted_calculator")
        assert isinstance(
            summarizer.weighted_calculator, WeightedCompletenessCalculator
        )

    def test_calculate_order_metrics_integration(self, temp_database_dir):
        """Test integration of weighted metrics in calculate_order_metrics."""
        summarizer = FullSummarizer(temp_database_dir)

        # Create mock counts file
        counts_content = "OG21\t1\nOG22\t1\nOG23\t0\n"
        counts_file = temp_database_dir / "test_counts.txt"
        counts_file.write_text(counts_content)

        # Test the integration
        completeness, duplication, weighted_comp, confidence = (
            summarizer.calculate_order_metrics(counts_file, "TestOrder")
        )

        # Validate return values
        assert isinstance(completeness, float)
        assert isinstance(duplication, float)
        assert isinstance(weighted_comp, float)
        assert isinstance(confidence, float)

        # Traditional completeness should be 2/3 * 100 = 66.67%
        assert abs(completeness - 66.67) < 0.1

        # Weighted completeness should be different due to conservation weighting
        assert weighted_comp > 0
        assert confidence > 0

    def test_error_handling_in_integration(self, temp_database_dir):
        """Test error handling when weighted calculation fails."""
        summarizer = FullSummarizer(temp_database_dir)

        # Mock the weighted calculator to raise an exception
        with patch.object(
            summarizer.weighted_calculator,
            "calculate_weighted_completeness",
            side_effect=Exception("Mock error"),
        ):

            counts_file = temp_database_dir / "test_counts.txt"
            counts_file.write_text("OG21\t1\n")

            completeness, duplication, weighted_comp, confidence = (
                summarizer.calculate_order_metrics(counts_file, "TestOrder")
            )

            # Should fallback gracefully
            assert weighted_comp == completeness  # Fallback to traditional
            assert confidence == 50.0  # Default fallback confidence


class TestPerformanceAndScaling:
    """Test performance characteristics and scaling behavior."""

    def test_large_marker_set_performance(self, temp_marker_stats_file):
        """Test performance with large marker sets."""
        calculator = WeightedCompletenessCalculator(
            marker_stats_path=temp_marker_stats_file,
            weight_strategy="information_theoretic",
        )

        # Generate large marker set
        large_marker_set = [f"marker_{i}" for i in range(1000)]
        marker_counts = {
            marker: 1 if i % 3 != 0 else 0 for i, marker in enumerate(large_marker_set)
        }

        import time

        start_time = time.time()

        weighted_comp, confidence, metrics = calculator.calculate_weighted_completeness(
            marker_counts=marker_counts,
            taxonomic_order="test_order",
            expected_markers=large_marker_set,
        )

        execution_time = time.time() - start_time

        # Should complete reasonably quickly (< 1 second for 1000 markers)
        assert execution_time < 1.0
        assert isinstance(weighted_comp, float)
        assert isinstance(confidence, float)

    def test_caching_effectiveness(self, temp_marker_stats_file):
        """Test that caching improves performance on repeated calls."""
        calculator = WeightedCompletenessCalculator(
            marker_stats_path=temp_marker_stats_file, cache_size=64
        )

        markers = tuple(f"marker_{i}" for i in range(100))

        import time

        # First call (should be slower due to computation)
        start_time = time.time()
        weights1 = calculator.get_marker_weights("test_order", markers)
        first_call_time = time.time() - start_time

        # Second call (should be faster due to caching)
        start_time = time.time()
        weights2 = calculator.get_marker_weights("test_order", markers)
        second_call_time = time.time() - start_time

        # Results should be identical
        assert weights1 == weights2

        # Second call should be significantly faster (at least 2x)
        assert second_call_time < first_call_time / 2


def test_factory_function():
    """Test the create_weighted_calculator factory function."""
    with tempfile.TemporaryDirectory() as temp_dir:
        db_path = Path(temp_dir)

        # Create minimal marker stats file
        marker_stats = pd.DataFrame(
            {
                "marker_name": ["test_marker"],
                "order_name": ["test_order"],
                "percent_genomes_with_marker": [50.0],
            }
        )
        marker_stats.to_csv(db_path / "marker_stats.tsv", sep="\t", index=False)

        calculator = create_weighted_calculator(db_path)

        assert isinstance(calculator, WeightedCompletenessCalculator)
        assert calculator.weight_strategy == "information_theoretic"
        assert calculator.enable_adaptive_scaling is True


if __name__ == "__main__":
    # Configure logging for tests
    logging.basicConfig(level=logging.DEBUG)

    # Run tests
    pytest.main([__file__, "-v"])
