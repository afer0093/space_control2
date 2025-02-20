extern crate nalgebra as na;
extern crate rand;
extern crate ndarray;
extern crate rand_distr;

use na::{Vector3, Vector4, Matrix3, Matrix4, DMatrix, Quaternion};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use ndarray::{Array1, s};
use std::fs::File;
use std::io::{Write, BufWriter};

// 定数
const T_MAX: f64 = 100.0; // シミュレーション時間 [s]
const RPM_TO_RAD: f64 = std::f64::consts::PI / 30.0; // RPM → rad/s

// 慣性テンソル
const I_X: f64 = 2.5;
const I_Y: f64 = 1.5;
const I_Z: f64 = 2.0;
const INERTIA: Matrix3<f64> = Matrix3::new(I_X, 0.0, 0.0, 0.0, I_Y, 0.0, 0.0, 0.0, I_Z);
const INERTIA_INV: Matrix3<f64> = Matrix3::new(1.0 / I_X, 0.0, 0.0, 0.0, 1.0 / I_Y, 0.0, 0.0, 0.0, 1.0 / I_Z);

fn angular_velocity_derivative(omega: &Vector3<f64>, torque: &Vector3<f64>) -> Vector3<f64> {
    INERTIA_INV * (torque - omega.cross(&(INERTIA * omega)))
}

fn quaternion_derivative(q: &Quaternion<f64>, omega: &Vector3<f64>) -> Quaternion<f64> {
    let w = q.w;
    let x = q.i;
    let y = q.j;
    let z = q.k;

    let omega_x = omega.x;
    let omega_y = omega.y;
    let omega_z = omega.z;

    // Ω行列 (4x4)
    let omega_mat = Matrix4::new(
        0.0, omega_z, -omega_y, omega_x,
        -omega_z, 0.0, omega_x, omega_y,
        omega_y, -omega_x, 0.0, omega_z,
        -omega_x, -omega_y, -omega_z, 0.0,
    );

    // クォータニオンベクトル (q_x, q_y, q_z, q_w)
    let q_vec = Vector4::new(x, y, z, w);

    // dq/dt = (1/2) * Ω * q
    let dq_vec = 0.5 * omega_mat * q_vec;

    // クォータニオンの時間微分を返す
    Quaternion::new(dq_vec[3], dq_vec[0], dq_vec[1], dq_vec[2])
}

fn rk4_step(
    omega: &Vector3<f64>,
    q: &Quaternion<f64>,
    torque: &Vector3<f64>,
    dt: f64,
) -> (Vector3<f64>, Quaternion<f64>) {
    let k1_omega = angular_velocity_derivative(omega, torque);
    let k1_q = quaternion_derivative(q, omega);

    let k2_omega = angular_velocity_derivative(&(omega + k1_omega * (dt / 2.0)), torque);
    let k2_q = quaternion_derivative(
        &Quaternion::new(
            q.w + k1_q.w * (dt / 2.0),
            q.i + k1_q.i * (dt / 2.0),
            q.j + k1_q.j * (dt / 2.0),
            q.k + k1_q.k * (dt / 2.0),
        ),
        &(omega + k1_omega * (dt / 2.0)),
    );

    let k3_omega = angular_velocity_derivative(&(omega + k2_omega * (dt / 2.0)), torque);
    let k3_q = quaternion_derivative(
        &Quaternion::new(
            q.w + k2_q.w * (dt / 2.0),
            q.i + k2_q.i * (dt / 2.0),
            q.j + k2_q.j * (dt / 2.0),
            q.k + k2_q.k * (dt / 2.0),
        ),
        &(omega + k2_omega * (dt / 2.0)),
    );

    let k4_omega = angular_velocity_derivative(&(omega + k3_omega * dt), torque);
    let k4_q = quaternion_derivative(
        &Quaternion::new(
            q.w + k3_q.w * dt,
            q.i + k3_q.i * dt,
            q.j + k3_q.j * dt,
            q.k + k3_q.k * dt,
        ),
        &(omega + k3_omega * dt),
    );

    let new_omega = omega
        + (k1_omega + k2_omega * 2.0 + k3_omega * 2.0 + k4_omega) * (dt / 6.0);
    let new_q = Quaternion::new(
        q.w + (k1_q.w + k2_q.w * 2.0 + k3_q.w * 2.0 + k4_q.w) * (dt / 6.0),
        q.i + (k1_q.i + k2_q.i * 2.0 + k3_q.i * 2.0 + k4_q.i) * (dt / 6.0),
        q.j + (k1_q.j + k2_q.j * 2.0 + k3_q.j * 2.0 + k4_q.j) * (dt / 6.0),
        q.k + (k1_q.k + k2_q.k * 2.0 + k3_q.k * 2.0 + k4_q.k) * (dt / 6.0),
    ).normalize();
    
    (new_omega, new_q)
}

struct SpaceControl {
    delta_t: f64,
    i_x: f64,
    i_y: f64,
    i_z: f64,
    q: Matrix3<f64>,
    r: Matrix3<f64>,
    sigma_w: f64,
    sigma_v: f64,
    num_of_div: usize,
    num_of_kf_update: usize,
    q_true_0: Quaternion<f64>,
    omega_true_0: Vector3<f64>,
    q_est_0: Quaternion<f64>,
    omega_est_0: Vector3<f64>,
    p_0: DMatrix<f64>,
}

impl SpaceControl {
    fn estimate_state(&self, x_: &Array1<f64>, y_: &Array1<f64>, pbar: &DMatrix<f64>) -> (Array1<f64>, DMatrix<f64>) {
        let q_ = y_.slice(s![..4]);
        let omega_ = y_.slice(s![4..]);

        let a = DMatrix::from_row_slice(7, 7, &[
            0.0,                omega_[2] / 2.0,    -omega_[1] / 2.0,   omega_[0] / 2.0,    q_[3] / 2.0,                                    -q_[2] / 2.0,                                   q_[1] / 2.0,
            -omega_[2] / 2.0,   0.0,                omega_[0] / 2.0,    omega_[1] / 2.0,    q_[2] / 2.0,                                    q_[3] / 2.0,                                    -q_[0] / 2.0,
            omega_[1] / 2.0,    -omega_[0] / 2.0,   0.0,                omega_[2] / 2.0,    -q_[1] / 2.0,                                   q_[0] / 2.0,                                    q_[3] / 2.0,
            -omega_[0] / 2.0,   -omega_[1] / 2.0,   -omega_[2] / 2.0,   0.0,                -q_[0] / 2.0,                                   -q_[1] / 2.0,                                   -q_[2] / 2.0,
            0.0,                0.0,                0.0,                0.0,                0.0,                                            -(self.i_z - self.i_y) * omega_[2] / self.i_x,  -(self.i_z - self.i_y) * omega_[1] / self.i_x,
            0.0,                0.0,                0.0,                0.0,                -(self.i_x - self.i_z) * omega_[2] / self.i_y,  0.0,                                            -(self.i_x - self.i_z) * omega_[0] / self.i_y,
            0.0,                0.0,                0.0,                0.0,                -(self.i_y - self.i_x) * omega_[1] / self.i_z,  -(self.i_y - self.i_x) * omega_[0] / self.i_z,  0.0
        ]);

        let b = DMatrix::from_row_slice(7, 3, &[
            0.0,            0.0,            0.0, 
            0.0,            0.0,            0.0,
            0.0,            0.0,            0.0,
            0.0,            0.0,            0.0,
            1.0 / self.i_x, 0.0,            0.0,
            0.0,            1.0 / self.i_y, 0.0,
            0.0,            0.0,            1.0 / self.i_z
        ]);

        let phi = a.clone().scale(self.delta_t).exp();
        let gamma = match a.clone().try_inverse() {
            Some(inv_a) => inv_a * (phi.clone() - DMatrix::identity(7, 7)) * b,
            None => {
                // 擬似逆行列を使用
                let pseudo_inv_a = a.clone().pseudo_inverse(1e-6).unwrap_or_else(|_| {
                    eprintln!("Matrix pseudo-inversion failed in estimate_state");
                    DMatrix::zeros(7, 7)
                });
                pseudo_inv_a * (phi.clone() - DMatrix::identity(7, 7)) * b
            }
        };

        let x_next_ = phi.clone() * DMatrix::from_column_slice(7, 1, x_.as_slice().unwrap()); // ノイズ平均0
        let pbar_next = phi.clone() * pbar * phi.transpose() + gamma.clone() * DMatrix::from_fn(3, 3, |i, j| self.q[(i, j)] as f64) * gamma.transpose();

        let mut x_next_array = Array1::zeros(7);
        x_next_array.slice_mut(s![..4]).assign(&Array1::from_vec(x_next_.rows(0, 4).iter().cloned().collect::<Vec<f64>>()));
        x_next_array.slice_mut(s![4..]).assign(&Array1::from_vec(x_next_.rows(4, 3).iter().cloned().collect::<Vec<f64>>()));

        (x_next_array, pbar_next)
    }

    fn get_observation(&self, y_: &Array1<f64>, random_column: usize, v: &Vector3<f64>) -> Vector3<f64> {
        let q_ = y_.slice(s![..4]);

        let obs = DMatrix::from_row_slice(3, 3, &[
            q_[0].powi(2) - q_[1].powi(2) - q_[2].powi(2) + q_[3].powi(2),  2.0 * (q_[0] * q_[1] + q_[2] * q_[3]),                              2.0 * (q_[0] * q_[2] - q_[1] * q_[3]),
            2.0 * (q_[0] * q_[1] - q_[2] * q_[3]),                          -q_[0].powi(2) + q_[1].powi(2) - q_[2].powi(2) + q_[3].powi(2),     2.0 * (q_[1] * q_[2] + q_[0] * q_[3]),
            2.0 * (q_[0] * q_[2] + q_[1] * q_[3]),                          2.0 * (q_[1] * q_[2] - q_[0] * q_[3]),                              -q_[0].powi(2) - q_[1].powi(2) + q_[2].powi(2) + q_[3].powi(2),
        ]);

        let obs_column = obs.column(random_column);
        let obs_vec = Vector3::new(obs_column[0], obs_column[1], obs_column[2]);

        obs_vec + v
    }

    fn extended_kalman_filter(&self, x_: &Array1<f64>, z: &Vector3<f64>, y_: &Array1<f64>, pbar: &DMatrix<f64>, random_column: usize) -> (Array1<f64>, DMatrix<f64>) {
        let q_ = y_.slice(s![..4]);

        let h_temp = match random_column {
            0 => DMatrix::from_row_slice(3, 7, &[
                2.0 * q_[0],    -2.0 * q_[1],   -2.0 * q_[2],   2.0 * q_[3],   0.0, 0.0, 0.0,
                2.0 * q_[1],    2.0 * q_[0],    -2.0 * q_[3],   -2.0 * q_[2],  0.0, 0.0, 0.0,
                2.0 * q_[2],    2.0 * q_[3],    2.0 * q_[0],    2.0 * q_[1],   0.0, 0.0, 0.0
            ]),
            1 => DMatrix::from_row_slice(3, 7, &[
                2.0 * q_[1],    2.0 * q_[0],    2.0 * q_[3],    2.0 * q_[2],    0.0, 0.0, 0.0,
                -2.0 * q_[0],   2.0 * q_[1],    -2.0 * q_[2],   2.0 * q_[3],    0.0, 0.0, 0.0,
                -2.0 * q_[3],   2.0 * q_[2],    2.0 * q_[1],    -2.0 * q_[0],   0.0, 0.0, 0.0
            ]),
            2 => DMatrix::from_row_slice(3, 7, &[
                2.0 * q_[2],    -2.0 * q_[3],   2.0 * q_[0],    -2.0 * q_[1],   0.0, 0.0, 0.0,
                2.0 * q_[3],    2.0 * q_[2],    2.0 * q_[1],    2.0 * q_[0],    0.0, 0.0, 0.0,
                -2.0 * q_[0],   -2.0 * q_[1],   2.0 * q_[2],    2.0 * q_[3],    0.0, 0.0, 0.0
            ]),
            _ => DMatrix::zeros(3, 7), // デフォルトケース
        };

        let h = h_temp;

        let tmp = (h.clone() * pbar * h.transpose() + DMatrix::from_fn(3, 3, |i, j| self.r[(i, j)] as f64)).try_inverse().unwrap_or_else(|| {
            println!("Matrix inversion failed in extended_kalman_filter, using pseudo-inverse");
            (h.clone() * pbar * h.transpose() + DMatrix::from_fn(3, 3, |i, j| self.r[(i, j)] as f64)).pseudo_inverse(1e-10).unwrap_or_else(|_| {
                DMatrix::zeros(3, 3)
            })
        });
        let k = pbar * h.transpose() * tmp.clone();
        let p = pbar - k.clone() * h.clone() * pbar;

        let x_vec = DMatrix::from_column_slice(7, 1, x_.as_slice().unwrap());

        let x_next_ = x_vec.clone() + k.clone() * (z - h.clone() * x_vec.clone());

        let mut x_next_array = Array1::zeros(7);
        x_next_array.slice_mut(s![..4]).assign(&Array1::from_vec(x_next_.rows(0, 4).iter().cloned().collect::<Vec<f64>>()));
        x_next_array.slice_mut(s![4..]).assign(&Array1::from_vec(x_next_.rows(4, 3).iter().cloned().collect::<Vec<f64>>()));

        (x_next_array, p)
    }

    fn normalize_quaternion(&self, y_: &Array1<f64>) -> Array1<f64> {
        let q_ = y_.slice(s![..4]).to_owned(); // ここでスライスを所有権を持つArray1に変換
        let omega_ = y_.slice(s![4..]).to_owned(); // ここでスライスを所有権を持つArray1に変換

        let q = (q_[0].powi(2) + q_[1].powi(2) + q_[2].powi(2) + q_[3].powi(2)).sqrt();

        let q_normalized_ = q_ / q;

        let mut result = Array1::zeros(7);
        result.slice_mut(s![..4]).assign(&q_normalized_);
        result.slice_mut(s![4..]).assign(&omega_);

        result
    }

    fn run_simulation(&self) {
        let omega_true_0_array = Array1::from_vec(self.omega_true_0.iter().cloned().collect::<Vec<f64>>());
        let omega_est_0_array = Array1::from_vec(self.omega_est_0.iter().cloned().collect::<Vec<f64>>());
        let mut y_true_ = vec![{
            let mut arr = Array1::zeros(7);
            arr.slice_mut(s![..4]).assign(&Array1::from_vec(vec![self.q_true_0.i, self.q_true_0.j, self.q_true_0.k, self.q_true_0.w]));
            arr.slice_mut(s![4..]).assign(&omega_true_0_array);
            arr
        }];
        let mut y_est_ = vec![{
            let mut arr = Array1::zeros(7);
            arr.slice_mut(s![..4]).assign(&Array1::from_vec(vec![self.q_est_0.i, self.q_est_0.j, self.q_est_0.k, self.q_est_0.w]));
            arr.slice_mut(s![4..]).assign(&omega_est_0_array);
            arr
        }];
        y_est_[0] = self.normalize_quaternion(&y_est_[0]);
        let mut x_ = vec![y_true_[0].clone() - y_est_[0].clone()];
        let mut pbar = self.p_0.clone();
        let mut p = vec![self.p_0.clone()];
        let mut x_updated = false;

        for i in 1..=self.num_of_div {
            // println!("Iteration: {}", i);
            let normal = Normal::new(0.0, self.sigma_w).unwrap();
            let w = Vector3::new(
                normal.sample(&mut rand::rng()),
                normal.sample(&mut rand::rng()),
                normal.sample(&mut rand::rng())
            );

            let y_true_last = y_true_.last().unwrap();
            let y_true_last_omega = Vector3::new(y_true_last[4], y_true_last[5], y_true_last[6]);
            let y_true_last_q = Quaternion::new(y_true_last[3], y_true_last[0], y_true_last[1], y_true_last[2]);

            let (y_true_next_omega, y_true_next_q) = rk4_step(&y_true_last_omega, &y_true_last_q, &w, self.delta_t);
            let mut y_true_next_ = Array1::zeros(7);
            y_true_next_.slice_mut(s![..4]).assign(&Array1::from_vec(vec![y_true_next_q.i, y_true_next_q.j, y_true_next_q.k, y_true_next_q.w]));
            y_true_next_.slice_mut(s![4..]).assign(&Array1::from_vec(vec![y_true_next_omega.x, y_true_next_omega.y, y_true_next_omega.z]));
            y_true_.push(y_true_next_);

            let y_est_last = y_est_.last().unwrap();
            let y_est_last_omega = Vector3::new(y_est_last[4], y_est_last[5], y_est_last[6]);
            let y_est_last_q = Quaternion::new(y_est_last[3], y_est_last[0], y_est_last[1], y_est_last[2]);

            let (y_est_next_omega, y_est_next_q) = rk4_step(&y_est_last_omega, &y_est_last_q, &Vector3::zeros(), self.delta_t); // 推定で外乱入れない
            let mut y_est_next_ = Array1::zeros(7);
            y_est_next_.slice_mut(s![..4]).assign(&Array1::from_vec(vec![y_est_next_q.i, y_est_next_q.j, y_est_next_q.k, y_est_next_q.w]));
            y_est_next_.slice_mut(s![4..]).assign(&Array1::from_vec(vec![y_est_next_omega.x, y_est_next_omega.y, y_est_next_omega.z]));
            y_est_.push(y_est_next_);

            if x_updated {
                let (x_next_, new_pbar) = self.estimate_state(&Array1::zeros(x_.last().unwrap().len()), &y_est_.last().unwrap(), &pbar);
                x_updated = false;
                x_.push(x_next_);
                pbar = new_pbar;
            } else {
                let (x_next_, new_pbar) = self.estimate_state(&x_.last().unwrap(), &y_est_.last().unwrap(), &pbar);
                x_.push(x_next_);
                pbar = new_pbar;
            }
            p.push(pbar.clone());

            if i % self.num_of_kf_update == 0 {
                let random_column = rand::rng().random_range(0..3);
                let normal = Normal::new(0.0, self.sigma_v).unwrap();
                let v = Vector3::new(
                    normal.sample(&mut rand::rng()),
                    normal.sample(&mut rand::rng()),
                    normal.sample(&mut rand::rng())
                );

                let y_true_last = y_true_.last().unwrap().clone();
                let y_est_last = y_est_.last().unwrap().clone();

                let z = self.get_observation(&y_true_last, random_column, &v) - self.get_observation(&y_est_last, random_column, &v);
                let (x_next_, p_next) = self.extended_kalman_filter(&x_.last().unwrap(), &z, &y_est_last, &pbar, random_column);
                y_est_.last_mut().unwrap().assign(&(y_est_last.clone() + x_next_.clone()));

                p.last_mut().unwrap().copy_from(&p_next);
                x_.last_mut().unwrap().assign(&x_next_);
                pbar = p_next;
                x_updated = true;
            }

            let y_true_last = y_true_.last().unwrap().clone();
            let y_est_last = y_est_.last().unwrap().clone();

            y_true_.last_mut().unwrap().assign(&self.normalize_quaternion(&y_true_last));
            y_est_.last_mut().unwrap().assign(&self.normalize_quaternion(&y_est_last));
        }

        // データをファイルに保存
        let file = File::create("simulation_data.dat").expect("Unable to create file");
        let mut writer = BufWriter::new(file);

        writeln!(writer, "# time true_qw true_qx true_qy true_qz true_omega_x true_omega_y true_omega_z est_qw est_qx est_qy est_qz est_omega_x est_omega_y est_omega_z error_qw error_qx error_qy error_qz error_omega_x error_omega_y error_omega_z sqrt_p_0 sqrt_p_1 sqrt_p_2 sqrt_p_3 sqrt_p_4 sqrt_p_5 sqrt_p_6").unwrap();
        for i in 0..y_true_.len() {
            let t = i as f64 * self.delta_t;
            let y_true = &y_true_[i];
            let y_est = &y_est_[i];
            let error = y_true - y_est;
            let sqrt_p: Vec<f64> = (0..7).map(|j| p[i][(j, j)].sqrt()).collect();
            writeln!(
                writer,
                "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                t,
                y_true[0], y_true[1], y_true[2], y_true[3], y_true[4], y_true[5], y_true[6],
                y_est[0], y_est[1], y_est[2], y_est[3], y_est[4], y_est[5], y_est[6],
                error[0], error[1], error[2], error[3], error[4], error[5], error[6],
                sqrt_p[0], sqrt_p[1], sqrt_p[2], sqrt_p[3], sqrt_p[4], sqrt_p[5], sqrt_p[6]
            ).unwrap();
        }
    }
}

fn main() {
    let omega_y_rpm = 17.0;
    let omega_y_rad = omega_y_rpm * RPM_TO_RAD;
    let omega_true_0 = Vector3::new(0.1, omega_y_rad + 0.1, 0.0);

    let mut rng = rand::rng();
    let q_est_0 = Quaternion::new(
        1.0 + rng.random::<f64>(),
        rng.random::<f64>(),
        rng.random::<f64>(),
        rng.random::<f64>(),
    ).normalize();
    let omega_est_0 = Vector3::new(
        0.1 + rng.random::<f64>(),
        omega_y_rad + 0.1 + rng.random::<f64>(),
        rng.random::<f64>(),
    );

    let sigma_w: f64 = 0.01;
    let sigma_v: f64 = 0.01;
    let delta_t: f64 = 0.01;

    let space_control = SpaceControl {
        delta_t: delta_t, // 時間更新間隔
        i_x: 2.5,
        i_y: 1.5,
        i_z: 2.0,
        q: Matrix3::identity() * sigma_w.powi(2),
        r: Matrix3::identity() * sigma_v.powi(2),
        sigma_w: sigma_w,
        sigma_v: sigma_v,
        num_of_div: (T_MAX / delta_t) as usize, // シミュレーション時間を時間更新間隔で割った値
        num_of_kf_update: (1.0 / delta_t) as usize, // 観測更新の間隔を時間更新間隔で割った値
        q_true_0: Quaternion::new(1.0, 0.0, 0.0, 0.0),
        omega_true_0: omega_true_0,
        q_est_0: q_est_0,
        omega_est_0: omega_est_0,
        p_0: DMatrix::identity(7, 7) * 100.0,
    };

    space_control.run_simulation();
}
