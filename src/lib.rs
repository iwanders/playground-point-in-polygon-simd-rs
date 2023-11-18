pub mod pnpoly;

pub mod print {
    use std::arch::x86_64::{__m256, __m256d, __m256i};
    #[allow(dead_code)]
    /// Print a vector of m128 type.
    pub fn pi(input: &__m256i) -> String {
        use std::arch::x86_64::_mm256_storeu_si256;
        let v: [u8; 32] = [0; 32];
        unsafe { _mm256_storeu_si256(v.as_ptr() as *mut _, *input) };
        format!("{:02X?}", v)
    }

    #[allow(dead_code)]
    /// Print a vector of m128d type.
    pub fn pf(input: &__m256) -> String {
        use std::arch::x86_64::_mm256_storeu_ps;
        let v: [f32; 8] = [0.0; 8];
        unsafe { _mm256_storeu_ps(v.as_ptr() as *mut _, *input) };
        format!("{:?}", v)
    }
    #[allow(dead_code)]
    /// Print a vector of m128d type.
    pub fn pd(input: &__m256d) -> String {
        use std::arch::x86_64::_mm256_storeu_pd;
        let v: [f64; 4] = [0.0; 4];
        unsafe { _mm256_storeu_pd(v.as_ptr() as *mut _, *input) };
        format!(
            "[{:?} ({:x?}), {:?} ({:x?}), {:?} ({:x?}), {:?} ({:x?})]",
            v[0],
            v[0].to_bits(),
            v[1],
            v[1].to_bits(),
            v[2],
            v[2].to_bits(),
            v[3],
            v[3].to_bits(),
        )
    }
}

const DO_PRINTS: bool = true;

#[allow(unused_macros)]
#[macro_export]
/// Helper print macro that can be enabled or disabled.
macro_rules! trace {
    () => (if crate::DO_PRINTS {println!("\n");});
    ($($arg:tt)*) => {
        if crate::DO_PRINTS {
            println!($($arg)*);
        }
    }
}
