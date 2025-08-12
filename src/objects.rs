use eqsolver::single_variable::FDNewton;
use crate::consts::{DOUBLE_PI, DUST_A_MIN, PC2CM, SIGMA, YEAR2SECOND};
use crate::utils::{search_interval, IntervalResult};

#[derive(Clone,Debug)]
pub struct UVSource {
    data_time: Vec<f64>,
    data_lumin: Vec<f64>,
    source_type: UVSourceType,
}

#[derive(Clone,Debug)]
pub enum UVSourceType {
    CONSTANT(f64),
    TDE,
    DRW,
    JUMP,
    USER,
}

impl UVSource {
    pub fn from_constant(val: f64) -> UVSource {
        UVSource {
            data_time: vec![0.0],
            data_lumin: vec![val],
            source_type: UVSourceType::CONSTANT(val),
        }
    }

    pub fn new(data_time: Vec<f64>, data_lumin: Vec<f64>) -> UVSource {
        UVSource {
            data_time,
            data_lumin,
            source_type: UVSourceType::USER,
        }
    }

    pub fn get_value(&self, t: f64) -> f64 {
        // 根据不同的 UVSourceType 进行匹配
        match self.source_type {
            UVSourceType::CONSTANT(val) => val,
            _ => self.get_value_general_type(t),
        }
    }

    fn get_value_general_type(&self, t: f64) -> f64 {
        let data_time = &self.data_time;
        let data_lumin = &self.data_lumin;
        let size = data_time.len();

        let interval = search_interval(data_time, t);
        let idx = match interval {
            // 成功寻找到插值区间
            IntervalResult::Found(idx) => idx,

            // 边界处理
            IntervalResult::BelowFirst => return data_lumin[0],
            IntervalResult::AboveLast => return data_lumin[size-1],
            IntervalResult::InsufficientEdges => return data_lumin[0],
        };

        // 获取插值区间
        let (t0, t1) = (data_time[idx], data_time[idx+1]);
        let (lumin0, lumin1) = (data_lumin[idx], data_lumin[idx+1]);

        // 线性插值
        let f = (t - t0) / (t1 - t0);
        lumin0 + f * (lumin1 - lumin0)
    }
}


#[derive(Clone,Debug)]
pub struct DustCube {
    // position
    r: f64,      // pc
    theta: f64,  // rad
    phi: f64,    // rad

    // cube size
    d_r: f64,
    d_theta: f64,
    d_phi: f64,

    // physical parameters
    n_dust: f64,       // number density of dust in cm-3
    n_gas: f64,        // number density of gas in cm-3
    a: f64,            // dust size in cm
    temperature: f64,  // K
    tau_uv: f64,       // dimensionless
    tau_ir: f64,       // dimensionless
}

impl DustCube {
    pub fn new(
        coord: (f64, f64, f64),
        cube_size: (f64, f64, f64),
        n_dust: f64,
        n_gas: f64,
        dust_size: f64,
        temperature: f64,
    ) -> Self {
        // TODO check parameters

        let mut cube = DustCube {
            r: coord.0, theta: coord.1, phi: coord.2,
            d_r: cube_size.0, d_theta: cube_size.1, d_phi: cube_size.2,
            n_dust,
            n_gas,
            a: dust_size,
            tau_ir: 0.0,
            tau_uv: 0.0,
            temperature,
        };

        cube.update_tau_uv();

        cube
    }

    pub fn build_cubes_in_cone(
        // the edges of radii
        edge_radius: Vec<f64>,
        orientation: (f64, f64),
        cone_size: (f64, f64),
        // physical parameters
        n_dust: f64, n_gas: f64, dust_size: f64, temperature: f64
    ) -> Vec<Self> {
        let num_cubes = edge_radius.len() - 1;
        let mut cubes = Vec::<Self>::with_capacity(num_cubes);

        for i in 1..num_cubes {
            let r1 = edge_radius[i-1];
            let r2 = edge_radius[i];
            let r = (r1 + r2) / 2.0;
            let d_r = r2 - r1;

            let coord = (r, orientation.0, orientation.1);
            let cube_size = (d_r, cone_size.0, cone_size.1);
            let cube = DustCube::new(
                coord,
                cube_size,
                n_dust,
                n_gas,
                dust_size,
                temperature
            );
            cubes.push(cube);
        }

        cubes
    }

    pub fn is_free(&self) -> bool {
        self.a <= 0.0
    }

    // the unit of flux is erg/s/cm2
    pub fn absorb_energy(&mut self, flux: f64) {
        if self.is_free() {
            return;
        }

        let teq = |temp: f64| self.thermal_equilibrium(flux, temp);

        // Starting guess is 1500
        let s = FDNewton::new(teq).solve(1500.0);

        let temp = s.ok();
        match temp {
            Some(temp) => {
                self.temperature = temp;
            },
            None => {
                self.temperature = 0.0;
            }
        }
    }

    // dt: yr
    pub fn destroy_particle(&mut self, dt: f64) {
        let v = self.da_dt();
        let da = v * (dt * YEAR2SECOND);
        let new_a = self.a + da;

        if new_a < DUST_A_MIN {
            self.a = 0.0;
            self.temperature = 0.0;
        } else {
            self.a = new_a;
        }

        self.update_tau_uv();
    }

    /// 根据尘埃颗粒大小更新 tau_uv
    /// 该光深大小是径向的
    fn update_tau_uv(&mut self) {
        if self.is_free() {
            self.tau_uv = 0.0;
        }

        let a = self.a;
        let n_dust = self.n_dust;
        let d_r = self.d_r;

        self.tau_uv = DOUBLE_PI*a*a*n_dust*d_r*PC2CM;
    }

    fn thermal_equilibrium(&self, flux: f64, temp: f64) -> f64 {
        let q_uv = 1.;

        let a_5 = self.a / 1E-5;
        let k = 0.1 * a_5 * (temp/2300.);

        let q_abs = k / (k + 1.);

        let left = flux * q_uv;
        let right = q_abs*4.0*SIGMA*temp*temp*temp*temp
            + 8.328650681E19 * (-7E4/temp).exp();

        left - right
    }

    fn da_dt(&self) -> f64 {
        let temp = self.temperature;
        let coef = (-7E4/temp).exp();
        -coef * 1E15 * 1E23_f64.powf(-1.0/3.0)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;
    use crate::consts::{DUST_TEMPERATURE_MIN, QUAT_PI};
    use super::*;

    #[test]
    fn test_thermal_equilibrium() {
        let lumin = 2E45;

        let coord = (1.0, 0.0, 0.0);
        let cube_size = (0.5, PI/6.0, PI/6.0);
        let n_dust = 5E-9;
        let n_gas = 0.0;
        let dust_size = 1E-5;

        let mut cube = DustCube::new(
            coord, cube_size, n_dust, n_gas, dust_size, DUST_TEMPERATURE_MIN
        );

        let r = coord.0; // pc
        let flux = lumin / PC2CM / PC2CM / QUAT_PI / r / r;
        cube.absorb_energy(flux);
    }
}