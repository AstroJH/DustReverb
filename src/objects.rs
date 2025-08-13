use std::f64::consts::{FRAC_PI_2, PI};
use eqsolver::single_variable::FDNewton;
use crate::consts::{DOUBLE_PI, DUST_A_MIN, PC2CM, QUAT_PI, SIGMA, YEAR2SECOND};
use crate::utils::{search_interval, IntervalResult};
use crate::errors::{Result, Errors};

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

#[derive(Clone,Debug)]
pub struct DustCone {
    // orientation
    theta: f64,
    phi: f64,
    d_theta: f64,
    d_phi: f64,

    // the edges of radii
    edge_radius: Vec<f64>,

    // data
    cubes: Vec<DustCube>,
}

#[derive(Clone,Debug)]
pub struct DustShell {
    // the edges of thetas
    edge_theta: Vec<f64>,

    // the edges of phis
    edge_phi: Vec<f64>,

    // the edges of radius
    edge_radius: Vec<f64>,

    // data
    cones: Vec<DustCone>,
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

        let mut cube = Self {
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
        edge_radius: &Vec<f64>,
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

    pub fn get_dust_size(&self) -> f64 {
        self.a
    }

    pub fn get_tau_uv(&self) -> f64 {
        self.tau_uv
    }

    pub fn get_tau_ir(&self) -> f64 {
        self.tau_ir
    }

    pub fn get_temperature(&self) -> f64 {
        self.temperature
    }

    pub fn get_n_dust(&self) -> f64 {
        self.n_dust
    }

    pub fn get_n_gas(&self) -> f64 {
        self.n_gas
    }

    pub fn get_coord(&self) -> (f64, f64, f64) {
        (self.r, self.theta, self.phi)
    }

    pub fn get_cube_size(&self) -> (f64, f64, f64) {
        (self.d_r, self.d_theta, self.d_phi)
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

    pub fn emit(&self) -> f64 {
        panic!()
    }
}

impl DustCone {
    pub fn new(
        edge_radius: Vec<f64>,
        orientation: (f64, f64),
        cone_size: (f64, f64),
        n_dust: f64,
        n_gas: f64,
        dust_size: f64,
        temperature: f64
    ) -> Self {
        let cubes = DustCube::build_cubes_in_cone(
            &edge_radius, orientation, cone_size, n_dust, n_gas, dust_size, temperature
        );

        Self {
            theta: orientation.0,
            phi: orientation.1,
            d_theta: cone_size.0,
            d_phi: cone_size.1,
            edge_radius,
            cubes
        }
    }

    pub fn num_cube(&self) -> usize {
        self.cubes.len()
    }

    pub fn get_cube_mut(&mut self, i: usize) -> &mut DustCube {
        &mut self.cubes[i]
    }

    pub fn get_cube(&self, i: usize) -> &DustCube {
        &self.cubes[i]
    }

    pub fn locate_cube(&self, r: f64) -> Option<&DustCube> {
        let interval = search_interval(&self.edge_radius, r);

        match interval {
            IntervalResult::Found(idx) => Some(&self.cubes[idx]),
            _ => None
        }
    }

    pub fn locate_cube_mut(&mut self, r: f64) -> Option<&mut DustCube> {
        let interval = search_interval(&self.edge_radius, r);

        match interval {
            IntervalResult::Found(idx) => Some(&mut self.cubes[idx]),
            _ => None
        }
    }

    pub fn cubes_iter(&self) -> impl Iterator<Item=&DustCube> {
        self.cubes.iter()
    }

    pub fn update_temp_from_lumin(&mut self, lumin: f64, min_temp: f64) {
        // calculate optical depth
        let mut cumulative_tau = 0.0;
        let tau_uv_values: Vec<f64> = self.cubes
            .iter()
            .map(|cube| {
                let current_tau = cumulative_tau;
                cumulative_tau += cube.tau_uv;
                current_tau
            })
            .collect();

        for (cube, &tau_uv)
        in self.cubes.iter_mut().zip(tau_uv_values.iter()) {

            // calculate flux
            let d_squared = cube.r * cube.r; // pc x pc
            let flux =
                (lumin / PC2CM / PC2CM) / QUAT_PI / d_squared * (-tau_uv).exp();

            // update temperature
            cube.absorb_energy(flux);
            cube.temperature = cube.temperature.max(min_temp);
        }
    }

    pub fn test_and_debug(&self) {
        for (i, cube) in self.cubes.iter().enumerate() {
            println!(
                "Cube{}\tT={:.2} K\ttau_uv={:.2}\ta={:.2E}",
                i + 1,
                cube.temperature,
                cube.tau_uv,
                cube.a
            );
        }
    }
}

impl DustShell {
    // 各向同性、均匀质地的球形网格
    pub fn build_simple_shell(
        edge_radius: Vec<f64>,
        n_dust: f64,
        n_gas: f64,
        dust_size: f64,
        temperature: f64,
    ) -> Self {
        let orientation = (FRAC_PI_2, PI);
        let cone_size = (PI, DOUBLE_PI);
        let edge_theta = vec![0.0, PI];
        let edge_phi = vec![0.0, DOUBLE_PI];

        let cone = DustCone::new(
            edge_radius.clone(),
            orientation,
            cone_size,
            n_dust,
            n_gas,
            dust_size,
            temperature
        );

        Self {
            edge_theta,
            edge_phi,
            edge_radius,
            cones: vec![cone],
        }
    }

    pub fn cones_iter_mut(&mut self) -> impl Iterator<Item=&mut DustCone> {
        self.cones.iter_mut()
    }

    pub fn cones_iter(&self) -> impl Iterator<Item=&DustCone> {
        self.cones.iter()
    }

    pub fn locate_cube(&self, r: f64, theta: f64, phi: f64) -> Option<&DustCube> {
        // search the cone
        // [(theta1, phi1), (theta1, phi2), ..., (theta2, phi1), ...]
        let idx_theta = match search_interval(&self.edge_theta, theta) {
            IntervalResult::Found(idx) => idx,
            _ => return None
        };

        let idx_phi = match search_interval(&self.edge_phi, phi) {
            IntervalResult::Found(idx) => idx,
            _ => return None
        };

        let num_phi = self.edge_phi.len() - 1;
        let pos = idx_theta * num_phi + idx_phi;
        let cone = &self.cones[pos];


        // search the cube
        let edge_radius = &cone.edge_radius;
        let idx_radius = match search_interval(edge_radius, r) {
            IntervalResult::Found(idx) => idx,
            _ => return None
        };

        Some(&cone.cubes[idx_radius])
    }

    pub fn update_temp_from_lumin(&mut self, lumin: f64, min_temp: f64) {
        for cone in self.cones.iter_mut() {
            cone.update_temp_from_lumin(lumin, min_temp);
        }
    }

    pub fn test_and_debug(&self, i_cone: usize) {
        let cone = &self.cones[i_cone];
        cone.test_and_debug();
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;
    use crate::consts::DUST_TEMPERATURE_MIN;
    use super::*;

    #[test]
    fn test_uv_source() {
        let constant = 2E45;
        let source = UVSource::from_constant(constant);
        assert_eq!(source.get_value(0.5), constant);
        assert_eq!(source.get_value(-2.5), constant);
        assert_eq!(source.get_value(1.5), constant);
    }

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