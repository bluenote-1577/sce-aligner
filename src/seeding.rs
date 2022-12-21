use crate::types::*;
use std::time::Instant;
use std::arch::x86_64::*;
//use bio::data_structures::interval_tree::IntervalTree;
//use fxhash::{hash, FxHashMap, FxHashSet};
use rust_lapper::{Interval, Lapper};
use smallvec::SmallVec;

#[inline]
pub fn position_max<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value0.cmp(value1))
        .map(|(idx, _)| idx)
}

#[inline(always)]
fn min_val<T: Ord>(slice: &[T]) -> Option<&T> {
    slice
        .iter()
        .max_by(|value0, value1| value1.cmp(value0)).map(|x| x)
}

#[inline(always)]
fn position_min<T: Ord>(slice: &[T]) -> Option<usize> {
    slice
        .iter()
        .enumerate()
        .max_by(|(_, value0), (_, value1)| value1.cmp(value0))
        .map(|(idx, _)| idx)
}

#[inline]
#[target_feature(enable = "avx2")]
pub unsafe fn mm_hash256(kmer: __m256i) -> __m256i {
    let mut key = kmer;
    let s1 = _mm256_slli_epi64(key, 21);
    key = _mm256_add_epi64(key, s1);
    key = _mm256_xor_si256(key, _mm256_cmpeq_epi64(key, key));

    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 24));
    let s2 = _mm256_slli_epi64(key, 3);
    let s3 = _mm256_slli_epi64(key, 8);

    key = _mm256_add_epi64(key, s2);
    key = _mm256_add_epi64(key, s3);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 14));
    let s4 = _mm256_slli_epi64(key, 2);
    let s5 = _mm256_slli_epi64(key, 4);
    key = _mm256_add_epi64(key, s4);
    key = _mm256_add_epi64(key, s5);
    key = _mm256_xor_si256(key, _mm256_srli_epi64(key, 28));

    let s6 = _mm256_slli_epi64(key, 31);
    key = _mm256_add_epi64(key, s6);

    return key;
}

pub fn os_seeds(
    string: &[u8],
    k : usize,
    c: usize,
    contig_index: ContigIndex,
    new_sketch: &mut Sketch,
) {
    if new_sketch.kmer_seeds_k.is_none() {
        new_sketch.kmer_seeds_k = Some(KmerSeeds::default());
    }
    let kmer_seeds_k = &mut new_sketch.kmer_seeds_k;
    if string.len() < 2 * k {
        return;
    }
    let num_bits = std::mem::size_of::<SeedBits>() * 8;
    let s = k - c + 1;
    let reverse_shift_dist = 2 * k - 2;
    let max_mask = SeedBits::MAX >> (num_bits - 2 * k);
    let max_mask_s = SeedBits::MAX >> (num_bits - 2 * s);
    let max_rev_mask = !(0 | (3 << (2 * k - 2)));
    //    let max_rev_mask_aa = !(0 | (31 << (5 * (max_k_aa - 1))));

    let window_size = k - s + 1;
    let t = window_size / 2;
    //        dbg!(&orf);
    let phase = 0;
    let mut rolling_kmer_f: SeedBits = 0;
    let mut rolling_kmer_r: SeedBits = 0;
    let mut circ_buffer = vec![0; window_size];
    let mut running_circ = 0;

    let mut min_val_buf = u128::MAX;
    let mut min_pos = usize::MAX;

    for i in 0..string.len() - k + 1 {
        let nuc_f = BYTE_TO_SEQ[string[i] as usize];
        rolling_kmer_f <<= 2;
        rolling_kmer_f |= nuc_f;
        rolling_kmer_f &= max_mask;
        let nuc_r = 3 - nuc_f;
        rolling_kmer_r >>= 2;
        rolling_kmer_r &= max_rev_mask;
        rolling_kmer_r |= nuc_r << reverse_shift_dist;
        rolling_kmer_r &= max_mask;
        
        if i >= s - 1 {
            let rolling_smer_f = rolling_kmer_f & max_mask_s;
            let rolling_smer_r = rolling_kmer_r & max_mask_s;
            let canon_smer = if rolling_smer_f < rolling_smer_r {
                rolling_smer_f
            } else {
                rolling_smer_r
            };
            circ_buffer[running_circ] = mm_hash128(canon_smer as SeedBits);
            running_circ += 1;

            if i >= k - 1 {
                if min_pos == usize::MAX || min_pos == running_circ-1{
                    min_pos = position_min(&circ_buffer).unwrap();
                    min_val_buf = *min_val(&circ_buffer).unwrap();
                }
                else{
                    if circ_buffer[running_circ-1] < min_val_buf{
                        min_pos = running_circ-1;
                        min_val_buf = circ_buffer[running_circ-1];
                    }
                }
                if min_pos == (running_circ + t) % window_size {
                    let canonical;
                    let canon_kmer;
                    if rolling_kmer_f < rolling_kmer_r {
                        canon_kmer = rolling_kmer_f;
                        canonical = true;
                    } else {
                        canon_kmer = rolling_kmer_r;
                        canonical = false;
                    };
                    let kmer_positions = kmer_seeds_k.as_mut().unwrap()
                        .entry(canon_kmer)
//                        .or_insert(SmallVec::<[SeedPosition; 5]>::new());
                        .or_insert(vec![]);

                    kmer_positions.push(SeedPosition {
                        pos: i as GnPosition,
                        canonical: canonical,
                        contig_index,
                        phase: phase as u8,
                    });
                }
            }
        }
        if running_circ == window_size {
            running_circ = 0;
        }
    }
}
