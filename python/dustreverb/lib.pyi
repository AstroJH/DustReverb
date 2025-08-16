from numpy.typing import NDArray

def simulate(
    times: NDArray,
    edge_radius: NDArray,
    edge_theta: NDArray,
    edge_phi: NDArray,
    n_dust: float,
    dust_size: float,
    source_times: NDArray,
    source_lumin: NDArray,
    output_shell: str,
    output_curve: str
):
    """ Launch Dust Reverb Engine.

    Execute dust reverberation simulation to model the time-dependent response of dust
    surrounding a UV/optical source (e.g., active galactic nucleus) to central radiation.
    
    This function employs a spherical coordinate grid to compute the thermal re-radiation
    from dust grains after absorbing central UV/optical emission. The simulation outputs
    detailed dust response per grid shell and the integrated light curve in WISE band (W1/W2).
    
    Parameters
    ==========
    times : NDArray[float]
        Simulation timesteps (yr).

    edge_radius : NDArray[float]
        Radial grid boundaries (centimeters). Defines spherical shells for radial grid.
        Array length should be N_r + 1 where N_r is number of radial bins.
        Example: [r_min, r1, r2, ..., r_max]
    
    edge_theta : NDArray[float]
        Polar angle grid boundaries (radians). Defines theta-direction bins.
        Range: 0 (North pole) to PI (South pole). Length = N_theta + 1.
    
    edge_phi : NDArray[float]
        Azimuthal angle grid boundaries (radians). Defines phi-direction bins.
        Range: 0 to 2PI. Length = N_phi + 1.
    
    n_dust : float
        Dust number density (cm-3). Assumed uniform throughout the grid.
        Future versions will support spatially varying density.
    
    dust_size : float
        Characteristic dust grain size (cm). Used to compute absorption cross-section.
        Typical value: ~ 1e-5 cm (i.e., 0.1 um).
    
    source_times : NDArray[float]
        Central source light curve timestamps (yr).
        Defines the time variability of UV/optical radiation.
    
    source_lumin : NDArray[float]
        Central source luminosity at corresponding source_times (erg/s).
        Must match source_times length.
        Represents the heating radiation absorbed by dust.
    
    output_shell : str
        Output file path for shell records data.
    
    output_curve : str
        Output file path for MIR light curve.
    """
    ...

def calc_structure_function(): ...

# def 