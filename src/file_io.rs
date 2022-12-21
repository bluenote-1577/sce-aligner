use crate::seeding;
use crate::types::*;
use fxhash::FxHashMap;
use log::*;
use needletail::parse_fastx_file;
use rand::seq::SliceRandom;
use rand::thread_rng;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::sync::Mutex;

pub fn get_ref_size(ref_file: &str) -> usize{
    let reader = parse_fastx_file(&ref_file);
    let mut length = 0;
    if !reader.is_ok() {
        warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
    } else {
        let mut reader = reader.unwrap();
        while let Some(record) = reader.next() {
            let record = record.expect(&format!("Invalid record for file {}", ref_file));
            let seq = record.seq();
            length += seq.len()
        }
    }
    return length
}

pub fn fastx_to_sketches(ref_files: &Vec<String>, c: usize, k: usize) -> Vec<Sketch> {
    let ref_sketches: Mutex<Vec<_>> = Mutex::new(vec![]);
    let mut index_vec = (0..ref_files.len()).collect::<Vec<usize>>();
    index_vec.shuffle(&mut thread_rng());
    index_vec.into_par_iter().for_each(|i| {
        let ref_file = &ref_files[i];
        let mut new_sketch = Sketch::new(usize::MAX, c, k, ref_file.to_string(), false);
        let reader = parse_fastx_file(&ref_file);
        if !reader.is_ok() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        } else {
            let mut j = 0;
            let mut is_valid = false;
            let mut reader = reader.unwrap();
            while let Some(record) = reader.next() {
                if record.is_ok() {
                    let record = record.expect(&format!("Invalid record for file {}", ref_file));
                    let contig = record.id();
                    let seq = record.seq();
                    new_sketch
                        .contigs
                        .push(String::from_utf8(contig.to_vec()).unwrap());
                    debug!("Sketching {}", new_sketch.contigs[new_sketch.contigs.len()-1]);
                    new_sketch.contig_lengths.push(seq.len() as GnPosition);
                    new_sketch.contig_seq.push(seq.to_vec());

                    new_sketch.total_sequence_length += seq.len();
                    seeding::os_seeds(&seq, k, c, j as u32, &mut new_sketch);
                    //new_sketch.contig_order = 0;
                    j += 1;
                    is_valid = true;
                } else {
                    warn!("File {} is not a valid fasta/fastq file", ref_file);
                    is_valid = false;
                    break;
                }
            }
            if is_valid {
                {
                    let mut locked = ref_sketches.lock().unwrap();
                    locked.push(new_sketch);
                }
            }
        }
    });
    let mut ref_sketches = ref_sketches.into_inner().unwrap();
    ref_sketches.sort();
    return ref_sketches;
}
pub fn fastx_to_multiple_sketch_rewrite(
    ref_files: &Vec<String>,
    c: usize,
    k: usize,
) -> Vec<Sketch> {
    let ref_sketches: Mutex<Vec<_>> = Mutex::new(vec![]);
    let mut index_vec = (0..ref_files.len()).collect::<Vec<usize>>();
    index_vec.shuffle(&mut thread_rng());
    index_vec.into_par_iter().for_each(|i| {
        let ref_file = &ref_files[i];
        let reader = parse_fastx_file(&ref_file);
        if !reader.is_ok() {
            warn!("{} is not a valid fasta/fastq file; skipping.", ref_file);
        } else {
            let mut j = 0;
            let mut reader = reader.unwrap();
            trace!("Sketching {} {}", ref_file, i);
            while let Some(record) = reader.next() {
                if record.is_ok() {
                    let record = record.expect(&format!("Invalid record for file {}", ref_file));
                    let contig = record.id();
                    let seq = record.seq();
                    let mut new_sketch = Sketch::new(usize::MAX, c, k, ref_file.to_string(), false);
                    new_sketch
                        .contigs
                        .push(String::from_utf8(contig.to_vec()).unwrap());
                    new_sketch.contig_lengths.push(seq.len() as GnPosition);
                    new_sketch.contig_seq.push(seq.to_vec());

                    new_sketch.total_sequence_length += seq.len();
                    seeding::os_seeds(&seq, k, c , 0 as u32, &mut new_sketch);
                    new_sketch.contig_order = j;
                    let mut locked = ref_sketches.lock().unwrap();
                    locked.push(new_sketch);
                    j += 1;
                } else {
                    warn!("File {} is not a valid fasta/fastq file", ref_file);
                    break;
                }
            }
        }
    });
    let mut ref_sketches = ref_sketches.into_inner().unwrap();
    ref_sketches.sort();
    return ref_sketches;
}
