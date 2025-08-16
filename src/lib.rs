use numpy::{PyReadonlyArrayDyn, PyArrayMethods};
use pyo3::{pymodule, types::PyModule, PyResult, Python, Bound};

use std::fs::File;
use std::io::{BufWriter, Write};
use std::f64::consts::{FRAC_2_PI, PI};
use rayon::prelude::*;

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

fn simulate(
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


fn play(
    times: Vec<f64>,
    edge_radius: Vec<f64>,
    n_dust: f64,
    dust_size: f64,
    source_times: Vec<f64>,
    source_lumin: Vec<f64>,
    output_shell: &str,
    output_curve: &str
) {
    let shell = DustShell::build_simple_shell(
        edge_radius, n_dust, 0.0, dust_size, DUST_TEMPERATURE_MIN
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


pub fn calc_esf(
    lightcurves: &Vec<Lightcurve>,
    tau_lo: &[f64],
    tau_hi: &[f64]
) -> (Vec<f64>, Vec<f64>, Vec<usize>) {
    let mut res_tau = Vec::<f64>::new();
    let mut res_sf = Vec::<f64>::new();
    let mut res_num = Vec::<usize>::new();

    for (lo, hi) in tau_lo.iter().zip(tau_hi.iter()) {
        let (sum_tau, sum_diff, sum_diff_err_square, n) =
            parallel_sum_all_mag_diffs(lightcurves, *lo, *hi);

        let tau = sum_tau / n as f64;
        let diff = sum_diff / n as f64;
        let diff_err_square = sum_diff_err_square / n as f64;

        let sf = (FRAC_2_PI * diff * diff - diff_err_square).sqrt();
        res_tau.push(tau);
        res_sf.push(sf);
        res_num.push(n);
    }

    (res_tau, res_sf, res_num)
}


fn parallel_sum_all_mag_diffs(
    lightcurves: &Vec<Lightcurve>,
    tau_lo: f64,
    tau_hi: f64
) -> (f64, f64, f64, usize) {
    // lightcurves.par_iter();
    let (sum_tau, sum_diff, sum_diff_err_square, n): (f64, f64, f64, usize) =
        lightcurves.par_iter().map(|lc| {

            let (tau, diff, diff_err_square) =
                lc.find_mag_differences(tau_lo, tau_hi);
            let sum_tau = tau.iter().sum::<f64>();
            let sum_diff = diff.iter().sum::<f64>();
            let sum_diff_err_square = diff_err_square.iter().sum::<f64>();
            let n = tau.len();

            (sum_tau, sum_diff, sum_diff_err_square, n)
        }).reduce(
            || (0.0, 0.0, 0.0, 0),
            |(a1, a2, a3, a4), (b1, b2, b3, b4)| (a1+b1, a2+b2, a3+b3, a4+b4)
        );

    (sum_tau, sum_diff, sum_diff_err_square, n)
}


#[pymodule]
#[pyo3(name="lib")]
fn dust_reverb<'py>(_py: Python<'py>, m: &Bound<'py, PyModule>) -> PyResult<()> {
    // wrapper of `play`
    #[pyfn(m)]
    #[pyo3(name = "play")]
    fn play_py<'py>(
        times: PyReadonlyArrayDyn<'py, f64>,
        edge_radius: PyReadonlyArrayDyn<'py, f64>,
        n_dust: f64,
        dust_size: f64,
        source_times: PyReadonlyArrayDyn<'py, f64>,
        source_lumin: PyReadonlyArrayDyn<'py, f64>,
        output_shell: &str,
        output_curve: &str
    ) {
        let times = times.to_vec().unwrap();
        let edge_radius = edge_radius.to_vec().unwrap();
        let source_times = source_times.to_vec().unwrap();
        let source_lumin = source_lumin.to_vec().unwrap();

        play(
            times,
            edge_radius,
            n_dust,
            dust_size,
            source_times,
            source_lumin,
            output_shell,
            output_curve
        );
    }

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

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_play() {

    }
}