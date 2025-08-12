// x1 >= x0
pub fn linear_interpolation_point(
    x0: f64, x1: f64, y0: f64, y1: f64,
    target_x: f64) -> f64 {

    let x = vec![x0, x1];
    let y = vec![y0, y1];

    linear_interpolation(&x, &y, target_x).unwrap()
}

// 不考虑具有相同的自变量数据点
pub fn linear_interpolation(
    x: &Vec<f64>, y: &Vec<f64>, target_x: f64
) -> Option<f64> {
    let size = x.len();

    if size == 0 || size != y.len() {
        return None;
    }

    if size == 1 {
        return Some(y[0]);
    }

    let x_min = x[0];
    let x_max = x[size - 1];

    // 边界处理
    if target_x <= x_min {
        return Some(y[0]);
    }

    if target_x >= x_max {
        return Some(y[size-1]);
    }

    // 二分查找目标区间
    let idx = x.partition_point(
        |&val| val < target_x
    );

    // 获取插值区间
    let (x0, x1) = (x[idx-1], x[idx]);
    let (y0, y1) = (y[idx-1], y[idx]);

    // 计算插值比例
    let t = (target_x - x0) / (x1 - x0);

    // 线性插值公式
    Some(y0 + t * (y1 - y0))
}


pub enum IntervalResult {
    BelowFirst,
    AboveLast,
    Found(usize),
    InsufficientEdges,
}

// [left, right)
pub fn search_interval(edges: &Vec<f64>, target: f64) -> IntervalResult {
    if edges.len() < 2 {
        return IntervalResult::InsufficientEdges;
    }

    if target < edges[0] {
        return IntervalResult::BelowFirst;
    }

    let last_index = edges.len() - 1;
    if target >= edges[last_index] {
        return IntervalResult::AboveLast;
    }

    let mut low = 0;
    let mut high = last_index;

    while low < high {
        let mid = low + ((high - low) >> 1);

        if target < edges[mid] {
            high = mid;
        } else {
            low = mid + 1;
        }
    }

    IntervalResult::Found(high - 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_interpolation_except() {
        // 用户输入空数组
        let x = vec![];
        let y = vec![];
        let result = linear_interpolation(&x, &y, 1.0);
        assert_eq!(result, None);

        // 用户输入不等长的数组
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let y = vec![2.0, 3.0];
        let result = linear_interpolation(&x, &y, 1.0);
        assert_eq!(result, None);
    }

    #[test]
    fn test_linear_interpolation_work() {
        // 恰好找到对应的数据点
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![10.0, 2.0, 0.0, -1.0, 5.0];
        let result_y = linear_interpolation(&x, &y, 3.0);
        assert!(result_y.is_some());
        assert_eq!(result_y.unwrap(), 0.0);

        // 正常情况，进行线性插值
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![10.0, 2.0, 0.0, -1.0, 5.0];
        let result_y = linear_interpolation(&x, &y, 2.5);
        assert!(result_y.is_some());
        assert_eq!(result_y.unwrap(), 1.0);

        // 输入的 target_x 在自变量范围之外
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![10.0, 2.0, 0.0, -1.0, 5.0];
        let result_y = linear_interpolation(&x, &y, 0.0);
        assert!(result_y.is_some());
        assert_eq!(result_y.unwrap(), 10.0);

        let result_y = linear_interpolation(&x, &y, 6.0);
        assert!(result_y.is_some());
        assert_eq!(result_y.unwrap(), 5.0);

        // 只有一个数据点，视为常值函数
        let x = vec![0.0];
        let y = vec![7.0];
        let result_y = linear_interpolation(&x, &y, 3.0);
        assert!(result_y.is_some());
        assert_eq!(result_y.unwrap(), 7.0);
    }

    #[test]
    fn test_linear_interpolation_point() {
        let result_y = linear_interpolation_point(
            1.0, 1.0, 2.0, 3.0, 3.0
        );

        println!("{}", result_y);
    }
}