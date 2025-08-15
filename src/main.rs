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

fn input_times() -> Vec<f64> {
    vec![
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
        0.55, 0.6, 0.65, 0.7, 0.75, 0.80, 0.85, 0.9, 0.95, 1.0,
        1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5,
        1.55, 1.6, 1.65, 1.7, 1.75, 1.80, 1.85, 1.9, 1.95, 2.0,
        2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5,
        2.55, 2.6, 2.65, 2.7, 2.75, 2.80, 2.85, 2.9, 2.95, 3.0,
        3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5,
        3.55, 3.6, 3.65, 3.7, 3.75, 3.80, 3.85, 3.9, 3.95, 4.0,
        4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 4.5,
        4.55, 4.6, 4.65, 4.7, 4.75, 4.80, 4.85, 4.9, 4.95, 5.0,
        5.05, 5.1, 5.15, 5.2, 5.25, 5.3, 5.35, 5.4, 5.45, 5.5,
        5.55, 5.6, 5.65, 5.7, 5.75, 5.80, 5.85, 5.9, 5.95, 6.0,
        6.05, 6.1, 6.15, 6.2, 6.25, 6.3, 6.35, 6.4, 6.45, 6.5,
        6.55, 6.6, 6.65, 6.7, 6.75, 6.80, 6.85, 6.9, 6.95, 7.0,
        7.05, 7.1, 7.15, 7.2, 7.25, 7.3, 7.35, 7.4, 7.45, 7.5,
        7.55, 7.6, 7.65, 7.7, 7.75, 7.80, 7.85, 7.9, 7.95, 8.0,
        8.05, 8.1, 8.15, 8.2, 8.25, 8.3, 8.35, 8.4, 8.45, 8.5,
        8.55, 8.6, 8.65, 8.7, 8.75, 8.80, 8.85, 8.9, 8.95, 9.0,
    ]
}

fn main() {
    // TODO 提取命令行参数
    // 1. 输出文件
    // 2.

    // TODO 获取用户输入，构造 DustReverOption 实例
    let times = input_times();
    let rs = vec![0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4];

    let shell = DustShell::build_simple_shell(
        rs, 5E-9, 0.0, 1E-5, DUST_TEMPERATURE_MIN
    );

    // TODO 初始化 DustReverbEngine 实例
    let mut engine = DustReverbEngine::new(
        times.clone(),
        UVSource::new(
            vec![
                0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,
                1.5, 1.8, 2.4, 3.6, 4.3, 4.6, 5.3, 6.0
            ],
            vec![
                2E45, 3E45, 1E45, 9E44, 2.5E45, 3.5E45, 2E45, 3.4E45,
                4E45, 5E45, 3.8E45, 4.2E45, 5.1E45, 4.4E45, 5.4E45, 4.1E45
            ]
        ),
        // UVSource::from_jump(2E45, 4E45, 1.5),
        shell,
        DUST_TEMPERATURE_MIN
    );

    // TODO 启动 DustReverbEngine
    engine.fire();

    // 将数据保存至磁盘
    // 尘埃辐射可转至 Python 进行计算
    engine.shutdown("data.csv").unwrap();


    let observer = (7.62E8, FRAC_2_PI, PI);

    let w1curve = engine.calc_lightcurve(WISE_W1_LAMBDA_CM, times.clone(), observer);
    let w2curve = engine.calc_lightcurve(WISE_W2_LAMBDA_CM, times.clone(), observer);

    let file = File::create("lcurve.csv").unwrap();
    let mut writer = BufWriter::new(file);
    writeln!(&mut writer, "time,w1flux,w2flux").unwrap();

    for i in 0..times.len() {
        let time = times[i];
        let w1flux = w1curve[i];
        let w2flux = w2curve[i];

        writeln!(&mut writer, "{},{},{}", time, w1flux, w2flux).unwrap();
    }

}
