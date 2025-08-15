use numpy::{PyReadonlyArrayDyn, PyArrayMethods};
use pyo3::{pymodule, types::PyModule, PyResult, Python, Bound};

use std::fs::File;
use std::io::{BufWriter, Write};
use std::f64::consts::{FRAC_2_PI, PI};

use crate::consts::{DUST_TEMPERATURE_MIN, WISE_W1_LAMBDA_CM, WISE_W2_LAMBDA_CM};
use crate::engine::DustReverbEngine;
use crate::objects::{DustShell, UVSource};

mod objects;
mod errors;
mod engine;
mod consts;
mod record;
mod utils;
mod options;


#[pymodule]
#[pyo3(name="lib")]
fn dust_reverb<'py>(_py: Python<'py>, m: &Bound<'py, PyModule>) -> PyResult<()> {
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

    Ok(())
}