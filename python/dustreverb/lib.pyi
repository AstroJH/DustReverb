from numpy.typing import NDArray

def play(
    times: NDArray,
    edge_radius: NDArray,
    n_dust: float,
    dust_size: float,
    source_times: NDArray,
    source_lumin: NDArray,
    output_shell: str,
    output_curve: str
): ...