from pathlib import Path
import csv
from typing import Optional
import numpy as np


def generate_samples(
    output_csv: Path,
    num_windows: int,
    umbrella_min_start: float,
    umbrella_min_end: float,
    sample_prefix_name: Optional[str] = "sample",
) -> None:
    """
    Generate a samples.csv with columns
    """
    rows = []

    umbrella_centers = np.linspace(umbrella_min_start, umbrella_min_end, num_windows)

    for i, center in enumerate(umbrella_centers):
        print(i, center)
        run_number = i + 1
        center_str = f"{center:.4f}"

        sample_name = f"{sample_prefix_name}_{run_number}"

        rows.append(
            {
                "sample_name": sample_name,
                "umbrella_center": center_str,
                "run_number": run_number,
            }
        )

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "sample_name",
                "umbrella_center",
                "run_number",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    # output CSV file
    output_csv = Path("config/samples.csv")

    sample_prefix_name = "window"

    # umbrella windows
    num_windows = 25
    umbrella_min_start = (
        -28  # each window defines the minimum of the umbrella potential Å
    )
    umbrella_min_end = 28  # each window defines the minimum of the umbrella potential Å
    # ======================================

    generate_samples(
        output_csv,
        num_windows,
        umbrella_min_start,
        umbrella_min_end,
        sample_prefix_name,
    )
