use crate::consts::{LIGHT_SPEED, PC2CM, QUAT_PI, YEAR2SECOND, CGS2JANSKY};
use crate::objects::{DustCone, DustCube, DustShell, UVSource};
use crate::record::DiscreteTimeMachine;
use crate::errors::Result;

pub struct DustReverbEngine {
    times: Vec<f64>,
    source: UVSource,
    shell: DustShell,
    time_machine: DiscreteTimeMachine,
}

impl DustReverbEngine {
    pub fn new(
        times: Vec<f64>,
        source: UVSource,
        mut shell: DustShell,
        temp_min: f64
    ) -> Self {
        // TODO check times

        let t0 = times[0];
        let lumin0 = source.get_value(t0);
        let mut time_machine = DiscreteTimeMachine::new();
        shell.update_temp_from_lumin(lumin0, temp_min);

        time_machine.record(t0, &shell);

        DustReverbEngine {
            times,
            source,
            shell,
            time_machine,
        }
    }

    pub fn fire(&mut self) {
        let times = &self.times;
        let source = &self.source;
        let shell = &mut self.shell;
        let time_machine = &mut self.time_machine;

        for i in 0..times.len()-1 {
            let t_curr = times[i];
            let t_next = times[i+1];
            // let d_t = t_next - t_curr;

            for cone in shell.cones_iter_mut() {
                Self::update_next_frame_of_cone(
                    source, time_machine, cone, t_curr, t_next
                );
            }

            time_machine.record(t_next, &shell);
        }
    }

    pub fn shutdown(&self, output: &str) -> Result<()> {
        self.time_machine.write(output)
    }

    fn update_next_frame_of_cone(
        source: &UVSource,
        time_machine: &DiscreteTimeMachine,
        cone: &mut DustCone,
        t_cur: f64, t_next: f64
    ) {
        // 计算每一个 cube 的累积推迟光深
        let num_cube = cone.num_cube();
        let mut tau_uv_tr  = Vec::<f64>::with_capacity(num_cube);

        for i in 0..num_cube {
            let cube = cone.get_cube(i);

            let mut tau = 0.0;
            for ii in 0..i {
                let inner_cube = cone.get_cube(ii);

                let d = cube.get_coord().0 - inner_cube.get_coord().0; // pc
                let tr = t_next - d*PC2CM / LIGHT_SPEED / YEAR2SECOND; // yr


                let inner_cube_coord = inner_cube.get_coord();
                let inner_cube_record =
                    time_machine.get_cube(
                        tr,
                        inner_cube_coord.0,
                        inner_cube_coord.1,
                        inner_cube_coord.2
                    );

                // 如果 inner_cube_record 是 None ...
                match inner_cube_record {
                    Some(inner_cube_record) => {
                        tau += inner_cube_record.get_tau_uv();
                    },
                    None => {
                        let inner_cube_first_record = time_machine.get_first_cube_record(
                            inner_cube_coord.0,
                            inner_cube_coord.1,
                            inner_cube_coord.2
                        ).unwrap();

                        tau += inner_cube_first_record.get_tau_uv();
                    }
                }
            }

            tau_uv_tr.push(tau);
        }

        // 计算 cone 中的每一个 cube
        for i in 0..num_cube {
            let cube = cone.get_cube_mut(i);

            Self::calc_frame_in_cube(source, cube, t_cur, t_next, tau_uv_tr[i]);
        }
    }

    // calculate the state of cube at next time
    fn calc_frame_in_cube(
        source: &UVSource, cube: &mut DustCube,
        t_cur: f64, t_next: f64, tau_uv: f64
    ) {
        let dt = t_next - t_cur; // yr
        let distance = cube.get_coord().0; // pc

        // 计算抵达 cube 的流量
        let tr = t_next - (distance*PC2CM/LIGHT_SPEED)/YEAR2SECOND;
        let lumin = source.get_value(tr); // erg/s
        let flux =
            (lumin / PC2CM / PC2CM) / distance / distance / QUAT_PI * (-tau_uv).exp();

        // 升华尘埃
        // 用当前时刻的温度更新下一时刻的尘埃半径 a
        // 光深与 a 相关，被自动更新
        cube.destroy_particle(dt);

        // 计算热平衡方程
        // 更新下一时刻的温度
        cube.absorb_energy(flux);
    }

    fn calc_separation(coord1: (f64, f64, f64), coord2: (f64, f64, f64)) -> f64 {
        let (r1, theta1, phi1) = coord1;
        let (r2, theta2, phi2) = coord2;

        let d2 = r1.powi(2) + r2.powi(2) - 2.0 * r1 * r2 * (
            theta1.sin() * theta2.sin() * (phi2 - phi1).cos()
            +
            theta1.cos() * theta2.cos()
        );

        d2.sqrt()
    }

    pub fn calc_lightcurve(&self, lam_cm: f64, times: Vec<f64>, observer: (f64, f64, f64)) -> Vec<f64> {
        let mut curve = Vec::<f64>::with_capacity(times.len());
        let time_machine = &self.time_machine;
        let shell = &self.shell;

        let (r_o, _, _) = observer;
        let t0 = r_o / LIGHT_SPEED * PC2CM / YEAR2SECOND;

        for t in times {
            let mut flux_t = 0.0;
            for cone in shell.cones_iter() {
                for cube in cone.cubes_iter() {
                    let coord = cube.get_coord();
                    let (r, theta, phi) = coord;
                    let distance = Self::calc_separation(coord, observer); // pc

                    let t_tr = t - distance / LIGHT_SPEED * PC2CM / YEAR2SECOND + t0;
                    
                    let cube_record = time_machine.get_cube(t_tr, r, theta, phi);
                    let cube_record = match cube_record {
                        Some(obj) => obj,
                        None => time_machine.get_first_cube_record(r, theta, phi).unwrap()
                    };

                    let lumin = cube_record.emit(lam_cm);
                    let flux = lumin / PC2CM / PC2CM * CGS2JANSKY / QUAT_PI / distance / distance;

                    flux_t += flux;
                }
            }

            curve.push(flux_t);
        }

        curve
    }
}