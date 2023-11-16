pub mod pnpoly;

pub mod print {
    use std::arch::x86_64::{__m128, __m128d, __m128i};
    #[allow(dead_code)]
    /// Print a vector of m128 type.
    pub fn pi(input: &__m128i) -> String {
        use std::arch::x86_64::_mm_storeu_si128;
        let v: [u8; 16] = [0; 16];
        unsafe { _mm_storeu_si128(v.as_ptr() as *mut _, *input) };
        format!("{:02X?}", v)
    }

    #[allow(dead_code)]
    /// Print a vector of m128d type.
    pub fn pf(input: &__m128) -> String {
        use std::arch::x86_64::_mm_storeu_ps;
        let v: [f32; 4] = [0.0; 4];
        unsafe { _mm_storeu_ps(v.as_ptr() as *mut _, *input) };
        format!("{:?}", v)
    }
    #[allow(dead_code)]
    /// Print a vector of m128d type.
    pub fn pd(input: &__m128d) -> String {
        use std::arch::x86_64::_mm_storeu_pd;
        let v: [f64; 2] = [0.0; 2];
        unsafe { _mm_storeu_pd(v.as_ptr() as *mut _, *input) };
        format!("{:?}", v)
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
