use std::{f64::consts::FRAC_2_PI, fs::File};
use std::sync::Arc;
use csv::Reader;
use rayon::prelude::*;
use rand::Rng;
use crate::errors::{Errors, Result};

#[derive(Debug, Clone)]
struct LightcurveData {
    mjd: Box<[f64]>,
    mag: Box<[f64]>,
    magerr: Box<[f64]>,
}

#[derive(Debug, Clone)]
pub struct Lightcurve {
    data: Arc<LightcurveData>, // Read-only
    pub redshift: f64
}

impl Lightcurve {
    pub fn new(
        mjd: Vec<f64>,
        mag: Vec<f64>,
        magerr: Vec<f64>,
        redshift: f64
    ) -> Self {
        assert_eq!(mjd.len(), mag.len(), "mjd and mag must have same length");
        assert_eq!(mjd.len(), magerr.len(), "mjd and magerr must have same length");

        Self {
            data: Arc::new(LightcurveData {
                mjd: mjd.into_boxed_slice(),
                mag: mag.into_boxed_slice(),
                magerr: magerr.into_boxed_slice()
            }),
            redshift,
        }
    }

    /// The length of light curve.
    pub fn len(&self) -> usize {
        self.data.mjd.len()
    }

    /// Load light curve data from file with csv format.
    pub fn load_from_csv(
        filename: &str, redshift: f64,
        mjd_name: &str, mag_name: &str, magerr_name: &str
    ) -> Result<Self> {
        let file = File::open(filename);

        if file.is_err() {
            return Err(Errors::FailedOpenFile);
        }

        let file = file.unwrap();
        let mut reader = Reader::from_reader(file);

        let headers = reader.headers().unwrap();
        let (mut i_mjd, mut i_mag, mut i_magerr) = (0usize, 0usize, 0usize);
        let mut flag = 0u8;

        for (i, header) in headers.iter().enumerate() {

            if header.eq(mjd_name) {
                i_mjd = i;
                flag += 1;
            } else if header.eq(mag_name) {
                i_mag = i;
                flag += 1;
            } else if header.eq(magerr_name) {
                i_magerr = i;
                flag += 1;
            }

            if flag >=3 {
                break;
            }
        }

        if flag < 3 {
            return Err(Errors::CsvHeaderError);
        }

        const CAPACITY: usize = 32;
        let mut mjd_vec = Vec::<f64>::with_capacity(CAPACITY);
        let mut mag_vec = Vec::<f64>::with_capacity(CAPACITY);
        let mut magerr_vec = Vec::<f64>::with_capacity(CAPACITY);

        for record in reader.records().into_iter() {
            let record = record.unwrap();
            let mjd = record.get(i_mjd).unwrap().parse::<f64>().unwrap();
            let mag = record.get(i_mag).unwrap().parse::<f64>().unwrap();
            let magerr = record.get(i_magerr).unwrap().parse::<f64>().unwrap();

            mjd_vec.push(mjd);
            mag_vec.push(mag);
            magerr_vec.push(magerr);
        }

        Ok(Lightcurve::new(mjd_vec, mag_vec, magerr_vec, redshift))
    }

    // TODO
    // Load light curve data from file with FITS format.
    // pub fn load_from_fits() {

    // }

    pub fn load_lightcurves(input: Vec<String>, redshift: Vec<f64>,
                            mjd_name: &str, mag_name: &str, magerr_name: &str) -> Vec<Self> {
        input.par_iter().zip(redshift.par_iter())
            .filter_map(|(filename, z)| {
                match Lightcurve::load_from_csv(filename, *z, mjd_name, mag_name, magerr_name) {
                    Ok(lightcurve) => Some(lightcurve),
                    Err(_) => None
                }
            })
            .collect()
    }

    /// 寻找时间间隔在 \[tau_lo, tau_hi\] 之间的点对，
    /// 返回这些点对的时间间隔、星等差以及对应的误差估计
    pub fn find_mag_differences(
        &self,
        tau_lo: f64,
        tau_hi: f64
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mjd = &self.data.mjd;
        let mag = &self.data.mag;
        let magerr = &self.data.magerr;
        let n = mjd.len();

        if n < 2 { return (Vec::new(), Vec::new(), Vec::new()) }

        // 计算相邻时间间隔
        let mut intervals = Vec::with_capacity(n-1);
        for i in 1..mjd.len() {
            intervals.push(
                (mjd[i] - mjd[i-1])/(1.0+self.redshift)
            );
        }

        let mut mag_diffs = Vec::<f64>::new();
        let mut mag_diff_errs = Vec::<f64>::new();
        let mut taus = Vec::<f64>::new();

        // lp := left pointer
        // rp := right pointer
        for lp in 0..(n-1) {
            let mut cumulative = intervals[lp];
            let mut rp = lp+1;

            while rp < n-1 && cumulative < tau_lo {
                cumulative += intervals[rp];
                rp += 1;
            }

            if cumulative < tau_lo {
                // 不再有时间间隔 >= tau_lo, stop!
                break;
            } else if cumulative > tau_hi {
                // 当前 right pointer 已不再有时间间隔位于 [tau_lo, tau_hi]
                // 下一轮尝试
                continue;
            }

            // 开始收集本轮满足要求的时间间隔点对
            while rp < n {
                let diff = (mag[lp] - mag[rp]).abs();
                let tau = mjd[rp] - mjd[lp];
                let diff_err = magerr[lp].powi(2) + magerr[rp].powi(2);

                mag_diffs.push(diff);
                taus.push(tau);
                mag_diff_errs.push(diff_err);

                if rp < n-1 {
                    cumulative += intervals[rp];
                    if cumulative > tau_hi {
                        break;
                    }
                }
                rp += 1;
            }
        }

        (taus, mag_diffs, mag_diff_errs)
    }

    /// 计算该光变曲线的结构函数
    pub fn calc_structure_function(&self, tau_lo: Vec<f64>, tau_hi: Vec<f64>)
    -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        
        let mut sf_tau = Vec::<f64>::with_capacity(16);
        let mut sf_val = Vec::<f64>::with_capacity(16);
        let mut sf_err = Vec::<f64>::with_capacity(16);
        
        for (lo, hi) in tau_lo.iter().zip(tau_hi.iter()) {
            let (tau, dmag, dmag_err)
            = self.find_mag_differences(*lo, *hi);
            
            // TODO
            let n = tau.len() as f64;
            let sum_tau: f64 = tau.iter().sum();
            sf_tau.push(sum_tau/n);
            // sf_val.push();
        }

        // ()

        (sf_tau, sf_val, sf_err)
    }

    /// 对光变曲线进行重采样
    pub fn resample(original: &[Self]) -> Vec<Self> {
        let n = original.len();
        let mut rng = rand::rng();

        (0..n)
            .map(|_| {
                let idx = rng.random_range(0..n);
                let lc = &original[idx];

                lc.clone()
            })
            .collect()
    }

    pub fn calc_esf(
        lightcurves: &Vec<Self>,
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
}

fn parallel_sum_all_mag_diffs(
    lightcurves: &Vec<Lightcurve>,
    tau_lo: f64,
    tau_hi: f64
) -> (f64, f64, f64, usize) {
    let (sum_tau, sum_diff, sum_diff_err_square, n) =
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



#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use crate::lightcurve::Lightcurve;
    use crate::errors::Result;

    fn load_a_sample() -> Result<Lightcurve> {
        let home = std::env::var("HOME").unwrap();
        let rtd = PathBuf::from(home)
            .join("Repository")
            .join("Paliya_2024_Seyfert1")
            .join("ztflc_repro");

        let filepath = rtd.join("10667-58163-0639.csv");

        let lcurve = Lightcurve::load_from_csv(
            filepath.to_str().unwrap(),
            0.3,
            "mjd",
            "mag",
            "magerr"
        );

        lcurve
    }

    #[test]
    fn test_ref() {
        let lcurve = load_a_sample().unwrap();
        let lc_clone = lcurve.clone();

        println!("{:p}", lcurve.data.mag);
        println!("{:p}", lc_clone.data.mag);
    }

    #[test]
    fn it_works() {
        let lcurve = load_a_sample().unwrap();
        println!("{:#?}", lcurve)
    }
}
