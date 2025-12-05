from snakemake.script import snakemake
import pandas as pd
from pathlib import Path
import json
from typing import Optional


def main(
    input_time_series_paths: list[Path],
    output_path: Path,
    umbrella_centers: list[float],
    spring_constants: list[float],
) -> None:
    """
    Create a WHAM meta file where each line has:
        <time_series_path> <umbrella_center> <spring_constant>
    """
    # Basic sanity checks
    n = len(input_time_series_paths)
    if len(umbrella_centers) != n:
        raise ValueError(
            f"umbrella_centers length ({len(umbrella_centers)}) "
            f"does not match number of time series paths ({n})"
        )

    # Allow a single spring constant for all windows, or one per window
    if len(spring_constants) == 1 and n > 1:
        springs = [float(spring_constants[0])] * n
    elif len(spring_constants) == n:
        springs = [float(k) for k in spring_constants]
    else:
        raise ValueError(
            f"spring_constants length ({len(spring_constants)}) must be 1 or {n}"
        )

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Write meta file
    with output_path.open("w") as f:
        for path, center, k in zip(input_time_series_paths, umbrella_centers, springs):
            # First column: each entry of abs_paths (as given by Snakemake)
            # Second column: umbrella center
            # Third column: spring constant
            f.write(f"{path} {center:g} {k:g}\n")


if __name__ == "__main__":
    main(
        snakemake.params.get("abs_paths"),  # type: ignore
        Path(snakemake.output.get("wham_meta")),  # type: ignore
        snakemake.params.get("umbrella_centers"),  # type: ignore
        snakemake.params.get("spring_constants"),  # type: ignore
    )
