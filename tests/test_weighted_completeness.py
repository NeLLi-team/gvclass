"""
Comprehensive tests for weighted completeness metrics.

Tests the ML-enhanced weighted completeness calculator with various scenarios
including edge cases, fallback behaviors, and performance characteristics.
"""

import logging
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from src.core.summarize import order_completeness
from src.core.summarize_full import FullSummarizer
from src.core.weighted_completeness import (
    WeightedCompletenessCalculator,
    create_weighted_calculator,
)


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

    from tests.conftest import skip_if_no_runtime_resources

    pytestmark = skip_if_no_runtime_resources()

    @pytest.fixture
    def temp_database_dir(self):
        """Create temporary database directory with required files."""
        from tests.conftest import stage_db_resources

        with tempfile.TemporaryDirectory() as temp_dir:
            db_path = Path(temp_dir)
            stage_db_resources(db_path, markers=False)

            # Create v1.5.0 nested layout for marker + order files.
            (db_path / "markers").mkdir()
            marker_stats = pd.DataFrame(
                {
                    "marker_name": ["OG21", "OG22", "OG23"],
                    "order_name": ["TestOrder"] * 3,
                    "percent_genomes_with_marker": [50.0, 75.0, 25.0],
                }
            )
            marker_stats.to_csv(db_path / "markers" / "stats.tsv", sep="\t", index=False)

            # Create mock order completeness file
            order_comp = pd.DataFrame(
                {
                    "Order": ["TestOrder"],
                    "Orthogroups": ["OG21, OG22, OG23"],
                    "Average_Percent": [66.67],
                    "Std_Percent": [5.0],
                }
            )
            order_comp.to_csv(db_path / "markers" / "order_completeness.tab", sep="\t", index=False)

            # Create mock labels file
            labels_content = "genome1\tNCLDV|Class|TestOrder|Family|Genus|Species\n"
            (db_path / "labels.tsv").write_text(labels_content)

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

        marker_counts = {"OG21": 1, "OG22": 1, "OG23": 0}

        # Test the integration
        metrics = order_completeness.calculate_order_metrics(
            summarizer.completeness_table,
            summarizer.order_stats_df,
            summarizer.weighted_calculator,
            summarizer.novelty_scorer,
            summarizer.completeness_mode,
            marker_counts,
            "TestOrder",
        )

        # Validate return values
        assert isinstance(metrics["order_completeness"], float)
        assert isinstance(metrics["order_dup"], float)
        assert isinstance(metrics["order_weighted_completeness"], float)
        assert isinstance(metrics["order_confidence_score"], float)

        # Traditional completeness should be 2/3 * 100 = 66.67%
        assert abs(metrics["order_completeness_raw"] - 66.67) < 0.1
        assert abs(metrics["order_completeness"] - 100.0) < 0.1

        # Weighted completeness should be different due to conservation weighting
        assert metrics["order_weighted_completeness"] > 0
        assert metrics["order_confidence_score"] > 0

    def test_failing_weighted_calculation_is_not_hidden(self, temp_database_dir):
        """A weighted calculation that fails raises instead of reporting a fallback."""
        summarizer = FullSummarizer(temp_database_dir)

        with (
            patch.object(
                summarizer.weighted_calculator,
                "calculate_weighted_completeness",
                side_effect=ValueError("Mock error"),
            ),
            pytest.raises(ValueError, match="Mock error"),
        ):
            order_completeness.calculate_order_metrics(
                summarizer.completeness_table,
                summarizer.order_stats_df,
                summarizer.weighted_calculator,
                summarizer.novelty_scorer,
                summarizer.completeness_mode,
                {"OG21": 1},
                "TestOrder",
            )


class TestPerformanceAndScaling:
    """Test performance characteristics and scaling behavior."""

    def test_large_marker_set_is_scored(self, temp_marker_stats_file):
        """A marker set far larger than the statistics table is scored."""
        calculator = WeightedCompletenessCalculator(
            marker_stats_path=temp_marker_stats_file,
            weight_strategy="information_theoretic",
        )
        large_marker_set = [f"marker_{i}" for i in range(1000)]
        marker_counts = {
            marker: 1 if i % 3 != 0 else 0 for i, marker in enumerate(large_marker_set)
        }

        weighted_comp, confidence, _metrics = (
            calculator.calculate_weighted_completeness(
                marker_counts=marker_counts,
                taxonomic_order="test_order",
                expected_markers=large_marker_set,
            )
        )

        assert weighted_comp == pytest.approx(100.0 * 666 / 1000)
        assert 0.0 <= confidence <= 100.0

    def test_marker_weights_are_computed_once_per_instance(
        self, temp_marker_stats_file
    ):
        calculator = WeightedCompletenessCalculator(
            marker_stats_path=temp_marker_stats_file
        )
        markers = tuple(f"marker_{i}" for i in range(100))

        with patch.object(
            calculator,
            "_compute_marker_weights",
            wraps=calculator._compute_marker_weights,
        ) as compute:
            first = calculator.get_marker_weights("test_order", markers)
            second = calculator.get_marker_weights("test_order", markers)

        assert second is first
        assert compute.call_count == 1

    def test_marker_stats_are_read_once_per_instance(self, temp_marker_stats_file):
        calculator = WeightedCompletenessCalculator(
            marker_stats_path=temp_marker_stats_file
        )

        with patch.object(pd, "read_csv", wraps=pd.read_csv) as read_csv:
            first = calculator.marker_stats_df
            second = calculator.marker_stats_df

        assert second is first
        assert read_csv.call_count == 1

    def test_instances_do_not_share_marker_stats_or_weights(self, tmp_path):
        paths = []
        for name, percent in (("a", 40.0), ("b", 90.0)):
            stats_path = tmp_path / f"{name}.tsv"
            pd.DataFrame(
                {
                    "marker_name": ["OG1"],
                    "order_name": ["OrderA"],
                    "percent_genomes_with_marker": [percent],
                }
            ).to_csv(stats_path, sep="\t", index=False)
            paths.append(stats_path)
        first, second = (
            WeightedCompletenessCalculator(path, enable_adaptive_scaling=False)
            for path in paths
        )

        # Interleaved access: a cache shared by the instances would return the
        # statistics of whichever instance read first.
        first_df = first.marker_stats_df
        second_df = second.marker_stats_df
        first_weights = first.get_marker_weights("OrderA", ("OG1",))
        second_weights = second.get_marker_weights("OrderA", ("OG1",))

        assert first.marker_stats_df is first_df
        assert second.marker_stats_df is second_df
        assert first_df["percent_genomes_with_marker"].tolist() == [40.0]
        assert second_df["percent_genomes_with_marker"].tolist() == [90.0]
        assert first_weights != second_weights


def test_factory_function():
    """Test the create_weighted_calculator factory function."""
    with tempfile.TemporaryDirectory() as temp_dir:
        db_path = Path(temp_dir)

        # Create minimal marker stats file (v1.5.0 nested layout).
        (db_path / "markers").mkdir()
        marker_stats = pd.DataFrame(
            {
                "marker_name": ["test_marker"],
                "order_name": ["test_order"],
                "percent_genomes_with_marker": [50.0],
            }
        )
        marker_stats.to_csv(db_path / "markers" / "stats.tsv", sep="\t", index=False)

        calculator = create_weighted_calculator(db_path)

        assert isinstance(calculator, WeightedCompletenessCalculator)
        assert calculator.weight_strategy == "information_theoretic"
        assert calculator.enable_adaptive_scaling is True


if __name__ == "__main__":
    # Configure logging for tests
    logging.basicConfig(level=logging.DEBUG)

    # Run tests
    pytest.main([__file__, "-v"])


def test_weighted_completeness_fallback_averages_duplicate_marker_rows(
    tmp_path: Path,
) -> None:
    """Fallback rows of one marker in several orders weigh like their mean."""
    header = "order_name\tmarker_name\tpercent_genomes_with_marker\n"
    duplicated = tmp_path / "duplicated.tsv"
    duplicated.write_text(
        header + "OrderA\tOG1\t80.0\nOrderB\tOG1\t90.0\nOrderA\tOG2\t70.0\n"
    )
    averaged = tmp_path / "averaged.tsv"
    averaged.write_text(header + "OrderA\tOG1\t85.0\nOrderA\tOG2\t70.0\n")
    markers = ("OG1", "OG2")

    from_duplicates = WeightedCompletenessCalculator(
        duplicated, enable_adaptive_scaling=False
    ).get_marker_weights("UnknownOrder", markers)
    from_means = WeightedCompletenessCalculator(
        averaged, enable_adaptive_scaling=False
    ).get_marker_weights("UnknownOrder", markers)

    assert from_duplicates == pytest.approx(from_means)
