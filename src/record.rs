use std::fs::File;
use std::io::{BufWriter, Write};
use crate::objects::{DustCube, DustShell};
use crate::utils::{linear_interpolation_point, search_interval, IntervalResult};
use crate::errors::Result;

// 记录每一帧下 DustShell 的状态
#[derive(Debug)]
pub struct DiscreteTimeMachine {
    pub times: Vec<f64>,
    pub records: Vec<DustShell>,
}

impl DiscreteTimeMachine {
    pub fn new() -> DiscreteTimeMachine {
        DiscreteTimeMachine {
            times: vec![],
            records: vec![],
        }
    }

    pub fn record(&mut self, time: f64, shell: &DustShell) {
        self.times.push(time);
        self.records.push(shell.clone());
    }

    pub fn get_first_cube_record(&self, r: f64, theta: f64, phi: f64) -> Option<DustCube> {
        let times = &self.times;
        let records = &self.records;
        let size = times.len();

        if size == 0 {
            return None
        }

        match records[0].locate_cube(r, theta, phi) {
            Some(cube) => Some(cube.clone()),
            None => None
        }
    }

    pub fn get_last_cube_record(&self, r: f64, theta: f64, phi: f64) -> Option<DustCube> {
        let times = &self.times;
        let records = &self.records;
        let size = times.len();

        if size == 0 {
            return None
        }

        match records[size-1].locate_cube(r, theta, phi) {
            Some(cube) => Some(cube.clone()),
            None => None
        }
    }

    /// 获取指定时间、指定坐标的 cube。
    ///
    /// **Note**: 既不能追溯过去，不能得知未来。
    pub fn get_cube(&self, time: f64, r: f64, theta: f64, phi: f64)
                    -> Option<DustCube> {
        let times = &self.times;
        let records = &self.records;
        let size = times.len();

        if size == 0 {
            return None
        }

        // if size == 1 { ... }

        // current time
        if time == times[size-1] {
            let cube = records[size-1].locate_cube(r, theta, phi);
            return match cube {
                Some(cube) => Some(cube.clone()),
                None => None
            }
        }

        let interval_time = search_interval(times, time);

        let (idx0, idx1) = match interval_time {
            IntervalResult::Found(idx) => (idx, idx+1),
            _ => return None
        };

        let t0 = self.times[idx0];
        let t1 = self.times[idx1];

        let sh0 = &records[idx0];
        let sh1 = &records[idx1];

        /* TODO
            这里可以优化。在获取 cube0 时已经知道了 cube 的物理索引位置，
            不需要再重复根据坐标查询 cube1 的索引了。
         */
        let cube0 = match sh0.locate_cube(r, theta, phi) {
            Some(cube) => cube,
            None => return None
        };

        let cube1 = match sh1.locate_cube(r, theta, phi) {
            Some(cube) => cube,
            None => return None
        };

        let dust_size = linear_interpolation_point(
            t0, t1, cube0.get_dust_size(), cube1.get_dust_size(), time
        );

        let temperature = linear_interpolation_point(
            t0, t1, cube0.get_temperature(), cube1.get_temperature(), time
        );

        let coord = cube0.get_coord();
        let cube_size = cube0.get_cube_size();

        let cube = DustCube::new(
            coord, cube_size,
            cube0.get_n_dust(), cube0.get_n_gas(),
            dust_size, temperature
        );

        Some(cube)
    }

    pub fn write(&self, output: &str) -> Result<()> {
        let file = File::create(output).unwrap();

        let mut writer = BufWriter::new(file);
        writeln!(&mut writer, "time,r,theta,phi,temp,a,tau_uv,flux").unwrap();

        for (time, record) in self.times.iter().zip(self.records.iter()) {
            for cone in record.cones_iter() {
                for cube in cone.cubes_iter() {
                    let coord = cube.get_coord();
                    writeln!(&mut writer, "{},{},{},{},{},{},{}",
                             time, coord.0, coord.1, coord.2,
                             cube.get_temperature(),
                             cube.get_dust_size(),
                             cube.get_tau_uv(),
                    ).unwrap();
                }
            }
        }

        writer.flush().unwrap();
        Ok(())
    }
}