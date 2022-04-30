//! Provides the [error](https://en.wikipedia.org/wiki/Error_function) and related functions
use crate::algebra::abstr::Real;
use crate::elementary::Trigonometry;
use crate::special::gamma;
use crate::special::gamma::Gamma;
use crate::algebra::abstr::Sign;

/// Provides [error, also called the Gauss error function](https://en.wikipedia.org/wiki/Error_function) related functions
pub trait Error
{
    /// Error function
    ///
    /// ```math
    /// \operatorname{erf}(z) = \frac{2}{\sqrt\pi}\int_0^z e^{-t^2}\,dt.
    /// ```
    ///
    /// [For more information](https://en.wikipedia.org/wiki/Error_function)
    ///
    /// # Arguments
    ///
    /// self:
    ///
    fn erf(self) -> Self;

    /// Complementary error function
    ///
    /// ```math
    /// \operatorname{erfc}(z) = 1 - \operatorname{erf}(z),
    /// ```
    /// # Arguments
    ///
    /// self
    fn erfc(self) -> Self;

    /// Inverse error function
    ///
    /// ```math
    /// \operatorname{erfinv}(x) = \operatorname{erf}^{-1}(x) \quad \text{for} \quad x \in (-1, 1)
    /// ```
    ///
    /// [For more information](https://en.wikipedia.org/wiki/Error_function#Inverse_functions)
    ///
    /// # Arguments
    ///
    /// -1.0 < x < 1.0
    ///
    /// # Panics
    ///
    /// if x < -1.0 or x > 1.0
    ///
    fn erfinv(self) -> Self;

    /// Inverse complementary error function
    ///
    /// ```math
    /// \operatorname{erfc}^{-1}(1-z) = \operatorname{erf}^{-1}(z)
    /// ```
    /// # Arguments
    ///
    /// 0 <= x <= 2.0
    ///
    /// # Panics
    ///
    /// if x < 0.0 or x > 2.0
    ///
    fn erfcinv(self) -> Self;
}


impl Error for f64
{
    fn erf(self) -> Self
    {
        if self == 0.0
        {
            return self;
        }
        self.sign() * (1.0f64 - gamma::gamma_u(0.5f64, self * self) / Self::pi().sqrt())
    }

    fn erfc(self) -> Self
    {
        1.0 - self.erf()
    }

    /// <https://en.wikipedia.org/wiki/Error_function#Inverse_functions>
    /// Using the rational approximations tabulated in:
    ///J. M. Blair, C. A. Edwards, and J. H. Johnson,
    /// "Rational Chebyshev approximations for the inverse of the error function",
    fn erfinv(self) -> Self
    {
        let factors_leq075_p: [f64; 7] = [0.16030_49558_44066_229311e2,
            -0.90784_95926_29603_26650e2,
            0.18644_91486_16209_87391e3,
            -0.16900_14273_46423_82420e3,
            0.65454_66284_79448_7048e2,
            -0.86421_30115_87247_794e1,
            0.17605_87821_39059_0];

        let factors_leq075_q: [f64; 7] = [0.14780_64707_15138_316110e2,
            -0.91374_16702_42603_13936e2,
            0.21015_79048_62053_17714e3,
            -0.22210_25412_18551_32366e3,
            0.10760_45391_60551_23830e3,
            -0.20601_07303_28265_443e2, 0.1e1];

        let factors_leg09375_p: [f64; 8] = [-0.15238_92634_40726_128e-1,
            0.34445_56924_13612_5216,
            -0.29344_39867_25424_78687e1,
            0.11763_50570_52178_27302e2,
            -0.22655_29282_31011_04193e2,
            0.19121_33439_65803_30163e2,
            -0.54789_27619_59831_8769e1,
            0.23751_66890_24448];

        let factors_leg09375_q: [f64; 8] = [-0.10846_51696_02059_954e-1,
            0.26106_28885_84307_8511,
            -0.24068_31810_43937_57995e1,
            0.10695_12997_33870_14469e2,
            -0.23716_71552_15965_81025e2,
            0.24640_15894_39172_84883e2,
            -0.10014_37634_97830_70835e2,
            0.1e1];

        let factors_p: [f64; 9] = [0.10501_31152_37334_38116e-3,
            0.10532_61131_42333_38164_25e-1,
            0.26987_80273_62432_83544_516,
            0.23268_69578_89196_90806_414e1,
            0.71678_54794_91079_96810_001e1,
            0.85475_61182_21678_27825_185e1,
            0.68738_08807_35438_39802_913e1,
            0.36270_02483_09587_08930_02e1,
            0.88606_27392_96515_46814_9];

        let factors_q: [f64; 10] = [0.10501_26668_70303_37690e-3,
            0.10532_86230_09333_27531_11e-1,
            0.27019_86237_37515_54845_553,
            0.23501_43639_79702_53259_123e1,
            0.76078_02878_58012_77064_351e1,
            0.11181_58610_40569_07827_3451e2,
            0.11948_78791_84353_96667_8438e2,
            0.81922_40974_72699_07893_913e1,
            0.40993_87907_63680_15361_45e1,
            0.1e1];

        let a: f64 = self.abs();
        if a >= 1.0
        {
            if self == 1.0
            {
                return f64::INFINITY;
            }
            else
            {
                if self == -1.0
                {
                    return f64::NEG_INFINITY;
                }
            }
            panic!("|self| has to be <= 1.0")
        }
        else
        {
            if a <= 0.75
            {
                let t: f64 = self * self - 0.5625;

                self * horner(t, &factors_leq075_p) / horner(t, &factors_leq075_q)
            } else {
                if a <= 0.9375
                {
                    let t: f64 = self * self - 0.87890625;

                    self * horner(t, &factors_leg09375_p) / horner(t, &factors_leg09375_q)
                } else {
                    let t: f64 = 1.0 / (-(1.0 - a).ln()).sqrt();

                    horner(t, &factors_p) / (copysign(t, self) * horner(t, &factors_q))
                }
            }
        }
    }

    fn erfcinv(self) -> Self
    {
        if (self < 0.0f64) || (self > 2.0f64)
        {
            panic!("self is > 2.0 or < 0.0");
        }
        (1.0 - self).erfinv()
    }
}

impl Error for f32
{
    fn erf(self) -> Self
    {
        if self == 0.0
        {
            return self;
        }
        self.sign() * (1.0f32 - gamma::gamma_u(0.5f32, self * self) / Self::pi().sqrt())
    }

    fn erfc(self) -> Self
    {
        1.0 - self.erf()
    }

    fn erfinv(self) -> Self
    {
        let factors_leq075_p: [f32; 3] = [-0.13095_99674_222e2,
            0.26785_22576_022,
            -0.92890_57365e1];

        let factors_leq075_q: [f32; 4] = [-0.12074_94262_97e2,
            0.30960_61452_9e2,
            -0.17149_97799_1e2,
            0.1e1];

        let factors_leg09375_p: [f32; 4] = [-0.12402_56522_12e0,
            0.10688_05957_4e1,
            -0.19594_55607_8e1,
            0.42305_81357e0];

        let factors_leg09375_q: [f32; 4] = [0.88276_97997e-1,
            0.89007_43359e0,
            -0.21757_03119_6e1,
            0.1e1];

        let factors_p: [f32; 6] = [0.15504_70003_116e0,
            0.13827_19649_631e1,
            0.69096_93488_87e0,
            -0.11280_81391_617e1,
            0.68054_42468_25e0,
            -0.16444_15679_1e0];

        let factors_q: [f32; 3] = [0.15502_48498_22e0,
            0.13852_28141_995e1,
            0.1e1];

        let a: f32 = self.abs();
        if a >= 1.0
        {
            if self == 1.0
            {
                return f32::INFINITY
            } else {
                if self == -1.0
                {
                    return f32::NEG_INFINITY
                }
            }
            panic!("|self| has to be <= 1.0")
        } else {
            if a <= 0.75
            {
                let t: f32 = self * self - 0.5625;

                self * horner(t, &factors_leq075_p) / horner(t, &factors_leq075_q)
            } else {
                if a <= 0.9375
                {
                    let t: f32 = self * self - 0.87890625;

                    self * horner(t, &factors_leg09375_p) / horner(t, &factors_leg09375_q)
                } else {
                    let t: f32 = 1.0 / (-(1.0 - a).ln()).sqrt();

                    horner(t, &factors_p) / (copysign(t, self) * horner(t, &factors_q))
                }
            }
        }
    }

    fn erfcinv(self) -> Self
    {
        if !(0.0f32..=2.0f32).contains(&self)
        {
            panic!("self is > 2.0 or < 0.0");
        }
        (1.0 - self).erfinv()
    }
}

/// Error Function
///
/// ```math
/// \operatorname{erf}(z) = \frac{2}{\sqrt\pi}\int_0^z e^{-t^2}\,dt.
/// ```
///
/// [For more information](https://en.wikipedia.org/wiki/Error_function)
///
/// # Example
///
/// ```
/// use mathru::special::error;
///
/// let x: f64 = 0.0;
/// let error: f64 = error::erf(x);
/// ```
pub fn erf<T>(x: T) -> T
    where T: Real + Error + Gamma
{
    x.erf()
}

/// Complementary error function
///
/// ```math
/// \operatorname{erfc}(z) = 1 - \operatorname{erf} z,
/// ```
/// # Arguments
///
/// x
///
/// # Example
///
/// ```
/// use mathru::special::error;
///
/// let x: f64 = 0.0;
/// let error: f64 = error::erfc(x);
/// ```
pub fn erfc<T>(x: T) -> T
    where T: Real + Error + Gamma
{
    x.erfc()
}

/// Inverse error function
///
/// ```math
/// \operatorname{erfinv}(x) = \operatorname{erf}^{-1}(x) \quad \text{for} \quad x \in (-1, 1)
/// ```
///
/// # Arguments
///
/// -1.0 <= x <= 1.0
///
/// # Panics
///
/// if x < -1.0 or x > 1.0
///
/// # Example
///
/// ```
/// use mathru::special::error;
///
/// let x: f64 = 0.0;
/// let error: f64 = error::erfinv(x);
/// ```
pub fn erfinv<T>(x: T) -> T
    where T: Real + Error
{
    x.erfinv()
}

fn horner<T>(x: T, c: &[T]) -> T
    where T: Real
{
    let last_idx: usize = c.len() -1;
    let mut f: T = c[last_idx] * x;

    for i in (1..last_idx).rev()
    {
        f = (f + c[i]) * x;
    }

    f + c[0]
}

fn copysign<T>(a: T, b: T) -> T
    where T: Real
{
    let a_neg: bool = a < T::zero();
    let b_neg: bool = b < T::zero();
    if (a_neg && !b_neg) || (!a_neg && b_neg)
    {
        return -a;
    }
    a
}


/// Inverse complementary error function
///
/// ```math
/// \operatorname{erfc}^{-1}(1-z) = \operatorname{erf}^{-1}(z)
/// ```
/// # Arguments
///
/// 0 <= x <= 2.0
///
/// # Panics
///
/// if x < 0.0 or x > 2.0
///
/// # Example
///
/// ```
/// use mathru::special::error;
///
/// let x: f64 = 1.0;
/// let error: f64 = error::erfcinv(x);
/// ```
pub fn erfcinv<T>(x: T) -> T
    where T: Real + Error
{
    x.erfcinv()
}
