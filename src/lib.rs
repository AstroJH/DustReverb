use numpy::{IntoPyArray, PyArrayDyn, PyArrayMethods, PyReadonlyArrayDyn};
use numpy::ndarray::Array1;
use pyo3::{pymodule, PyResult, Python, Bound};
use pyo3::types::PyModule;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::f64::consts::{FRAC_2_PI, PI};


use crate::consts::{DUST_TEMPERATURE_MIN, WISE_W1_LAMBDA_CM, WISE_W2_LAMBDA_CM};
use crate::engine::DustReverbEngine;
use crate::lightcurve::Lightcurve;
use crate::objects::{DustShell, UVSource};

mod objects;
mod errors;
mod engine;
mod consts;
mod record;
mod utils;
mod options;
mod lightcurve;

/// 模拟尘埃升华
/// 
/// 
pub fn simulate(
    times: Vec<f64>,
    edge_radius: Vec<f64>,
    edge_theta: Vec<f64>,
    edge_phi: Vec<f64>,

    n_dust: f64,
    dust_size: f64,
    source_times: Vec<f64>,
    source_lumin: Vec<f64>,
    output_shell: &str,
    output_curve: &str
) {

    let shell = DustShell::from_dust(
        edge_radius,
        edge_theta,
        edge_phi,
        n_dust,
        0.0,
        dust_size,
        DUST_TEMPERATURE_MIN
    );

    let mut engine = DustReverbEngine::new(
        times.clone(),
        UVSource::new(source_times, source_lumin),
        shell,
        DUST_TEMPERATURE_MIN
    );

    engine.fire();

    engine.shutdown(output_shell).unwrap();

    let observer = (7.62E8, FRAC_2_PI, PI);

    let w1curve = engine.calc_lightcurve(WISE_W1_LAMBDA_CM, times.clone(), observer);
    let w2curve = engine.calc_lightcurve(WISE_W2_LAMBDA_CM, times.clone(), observer);

    let file = File::create(output_curve).unwrap();
    let mut writer = BufWriter::new(file);
    writeln!(&mut writer, "time,w1flux,w2flux").unwrap();

    for i in 0..times.len() {
        let time = times[i];
        let w1flux = w1curve[i];
        let w2flux = w2curve[i];

        writeln!(&mut writer, "{},{},{}", time, w1flux, w2flux).unwrap();
    }

}

pub fn calc_sf(
    mjd: Vec<f64>,
    mag: Vec<f64>,
    magerr: Vec<f64>,
    redshift: f64,
    tau_lo: Vec<f64>,
    tau_hi: Vec<f64>
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let lcurve = Lightcurve::new(mjd, mag, magerr, redshift);

    lcurve.calc_structure_function(tau_lo, tau_hi)
}

pub fn calc_esf(
    input: Vec<String>, redshifts: Vec<f64>,
    tau_lo: Vec<f64>, tau_hi: Vec<f64>,
    mjd_name: String, mag_name: String, magerr_name: String
) -> (Vec<f64>, Vec<f64>, Vec<usize>) {

    let lightcurves = Lightcurve::load_lightcurves(input, redshifts, &mjd_name, &mag_name, &magerr_name);
    println!("load {} light curves", lightcurves.len());
    let (tau, sf, num) = Lightcurve::calc_esf(&lightcurves, &tau_lo, &tau_hi);
    
    (tau, sf, num)
}

#[pymodule]
#[pyo3(name="lib")]
fn dust_reverb<'py>(_py: Python<'py>, m: &Bound<'py, PyModule>) -> PyResult<()> {
    #[pyfn(m)]
    #[pyo3(name = "simulate")]
    fn simulate_py<'py>(
        times: PyReadonlyArrayDyn<'py, f64>,
        edge_radius: PyReadonlyArrayDyn<'py, f64>,
        edge_theta: PyReadonlyArrayDyn<'py, f64>,
        edge_phi: PyReadonlyArrayDyn<'py, f64>,

        n_dust: f64,
        dust_size: f64,
        source_times: PyReadonlyArrayDyn<'py, f64>,
        source_lumin: PyReadonlyArrayDyn<'py, f64>,
        output_shell: &str,
        output_curve: &str
    ) {
        let times = times.to_vec().unwrap();
        let edge_radius = edge_radius.to_vec().unwrap();
        let edge_theta = edge_theta.to_vec().unwrap();
        let edge_phi = edge_phi.to_vec().unwrap();
        let source_times = source_times.to_vec().unwrap();
        let source_lumin = source_lumin.to_vec().unwrap();

        simulate(times, edge_radius, edge_theta, edge_phi, n_dust, dust_size, source_times, source_lumin, output_shell, output_curve);
    }

    #[pyfn(m)]
    #[pyo3(name = "calc_structure_function")]
    fn calc_sf_py<'py>(
        py: Python<'py>,
        mjd: PyReadonlyArrayDyn<'py, f64>,
        mag: PyReadonlyArrayDyn<'py, f64>,
        magerr: PyReadonlyArrayDyn<'py, f64>,
        redshift: f64,
        tau_lo: PyReadonlyArrayDyn<'py, f64>,
        tau_hi: PyReadonlyArrayDyn<'py, f64>
    ) -> (Bound<'py, PyArrayDyn<f64>>, Bound<'py, PyArrayDyn<f64>>, Bound<'py, PyArrayDyn<f64>>) {
        let mjd = mjd.to_vec().unwrap();
        let mag = mag.to_vec().unwrap();
        let magerr = magerr.to_vec().unwrap();
        let tau_lo = tau_lo.to_vec().unwrap();
        let tau_hi = tau_hi.to_vec().unwrap();

        let (tau, sf, err)
            = calc_sf(mjd, mag, magerr, redshift, tau_lo, tau_hi);
        
        (
            Array1::from_vec(tau).into_dyn().into_pyarray(py),
            Array1::from_vec(sf).into_dyn().into_pyarray(py),
            Array1::from_vec(err).into_dyn().into_pyarray(py)
        )
    }

    #[pyfn(m)]
    #[pyo3(name = "calc_ensemble_structure_function")]
    fn calc_esf_py<'py>(
        py: Python<'py>,
        input: Vec<String>,
        redshifts: PyReadonlyArrayDyn<'py, f64>,
        tau_lo: PyReadonlyArrayDyn<'py, f64>,
        tau_hi: PyReadonlyArrayDyn<'py, f64>,
        mjd_name: String,
        mag_name: String,
        magerr_name: String
    ) -> (Bound<'py, PyArrayDyn<f64>>, Bound<'py, PyArrayDyn<f64>>, Bound<'py, PyArrayDyn<usize>>) {

        let redshifts = redshifts.to_vec().unwrap();
        let tau_lo = tau_lo.to_vec().unwrap();
        let tau_hi = tau_hi.to_vec().unwrap();

        let (tau, sf, num) = calc_esf(input, redshifts, tau_lo, tau_hi, mjd_name, mag_name, magerr_name);

        (
            Array1::from_vec(tau).into_dyn().into_pyarray(py),
            Array1::from_vec(sf).into_dyn().into_pyarray(py),
            Array1::from_vec(num).into_dyn().into_pyarray(py)
        )
    }

    Ok(())
}

#[cfg(test)]
mod tests {
}