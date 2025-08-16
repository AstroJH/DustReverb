use std::f64::consts::PI;

pub const PC2CM: f64 = 3.086E18;
pub const YEAR2SECOND: f64 = 365.24 * 24.0 * 3600.0;
pub const CGS2JANSKY: f64 = 1E23;

pub const K_B: f64 = 1.380649e-16; // cgs erg/K
pub const H: f64 = 6.62607015e-27; // csg erg s
pub const SIGMA: f64 = 5.670373E-5; // cgs
pub const LIGHT_SPEED: f64 = 29979245800.; // cgs

pub const DOUBLE_PI: f64 = 2.0 * PI;
pub const TRIPLE_PI: f64 = 3.0 * PI;
pub const QUAT_PI: f64 = 4.0 * PI;

pub const DUST_A_MIN: f64 = 1E-6;
pub const DUST_TEMPERATURE_MIN: f64 = 150.0;

pub const WISE_W1_LAMBDA_CM: f64 = 0.00034;
pub const WISE_W2_LAMBDA_CM: f64 = 0.00046;
pub const FLUX_ZERO_W1: f64 = 306.681; // Jy
pub const FLUX_ZERO_W2: f64 = 170.663; // Jy