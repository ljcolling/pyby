fn welfords(
    (i, count, mean, squared_distances): &mut (usize, f64, f64, f64),
    new_value: &f64,
) -> Option<(usize, i32, f64, f64, f64)> {
    // Use Welford's method to work out the running mean and sd
    *count += 1f64;
    let delta = new_value - *mean;
    *mean += delta / *count;
    let delta2 = new_value - *mean;
    *squared_distances += delta * delta2;
    let sd = (*squared_distances / (*count - 1f64)).sqrt();

    Some((
        *i,
        *count as i32,     // item
        *new_value as f64, // current value
        *mean,             // running mean
        sd,                // running sd
                           // bf,
    ))
}
