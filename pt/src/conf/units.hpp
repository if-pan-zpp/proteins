#pragma once

inline long double operator""_AA(long double angstrem) {
    return angstrem / 5.0L;
}

inline long double operator""_dAA2(long double angstrem) {
    return angstrem * 5.0L * 5.0L;
}

inline long double operator""_dAA4(long double angstrem) {
    return angstrem * 5.0L * 5.0L * 5.0L * 5.0L;
}