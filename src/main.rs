use bio::alignment::pairwise::*;
use clap::{App, Arg, SubCommand};
use fxhash::FxHashSet;
use fxhash::FxHashMap;
use log::LevelFilter;
use rayon::prelude::*;
use sce_aligner::avl_tree::*;
use sce_aligner::file_io;
use sce_aligner::seeding;
use sce_aligner::types::*;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::mem;
use std::sync::Mutex;
use std::time::Instant;
fn main() {
    simple_logging::log_to_stderr(LevelFilter::Debug);
    let matches = App::new("basic seed chain aligner")
        .version("0.1")
        .about("Basic seed chain extend aligner on simulated sequences.")
        .arg(
            Arg::with_name("reference")
                .index(1)
                .required(true)
                .takes_value(true)
                .help("reference fasta file."),
        )
        .arg(
            Arg::with_name("query")
                .index(2)
                .required(true)
                .takes_value(true)
                .help("query fastx file."),
        )
        .arg(Arg::with_name("sketch").short('s').help(
            "Use sketching. Default density is 1/(k - 6) where k = C log n (default no sketching)",
        ))
        .arg(
            Arg::with_name("minimizer")
                .short('m')
                .help("Use minimizers instead of open syncmers"),
        )
        .arg(
            Arg::with_name("wfa")
                .long("wfa")
                .help("Use wavefront aligner instead of standard DP for extension")
                .hidden(true),
        )
        .arg(Arg::with_name("debug").long("debug").help("Debug mode"))
        .arg(
            Arg::with_name("threads")
                .short('t')
                .takes_value(true)
                .help("Number of threads (default 20)"),
        )
        .get_matches();
    let print_all_debug = false;
    let use_minimizers = matches.is_present("minimizer");
    let use_wfa = matches.is_present("wfa");
    let sketch = matches.is_present("sketch");
    let threads = matches
        .value_of("threads")
        .unwrap_or("20")
        .parse::<usize>()
        .unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();
    let n;
    let debug = matches.is_present("debug");
    let reference = matches.value_of("reference").unwrap().to_string();
    let mut ref_sketch = Sketch::default();
    if reference.contains(".sketch") {
        let reader = BufReader::new(File::open(&reference).expect(&reference));
        let res: Result<Sketch, _> = bincode::deserialize_from(reader);
        ref_sketch = res.unwrap();
        n = ref_sketch
            .contig_lengths
            .iter()
            .map(|x| (*x as usize))
            .sum();
    } else {
        n = file_io::get_ref_size(matches.value_of("reference").unwrap());
    }
    let theta = 0.05;
    let alpha = -((1.0 - theta) as f64).log(4.0);
    let k = ((n as f64).log(4.) * 2. / (1. - 2. * alpha)).round() as usize;
    let s;
    let sketch = true;
    if sketch {
        s = k - 6;
    } else {
        s = 1;
    }
    let w = 10.;
    let c = k - s + 1;
    let window = 2 * c - 1;
    let t = (k - s + 1) / 2 + 1;
    let zeta_inter = 2. / (1. - 2. * alpha) * 50. / 8.
        * (n as f64).log(4.)
        * (n as f64).ln()
        * ((1. - theta) as f64).powi(-(k as i32));
    let zeta_sketch_inter =
        zeta_inter * 4. * 3. / 2. + 8. * zeta_inter / ((1. - theta) as f64).powi(-(k as i32));
    let zeta_inter = 1. / (6. * zeta_inter);
    //            dbg!(6. * 1. / zeta);
    //            let zeta = 0.00;
    let zeta_sketch_inter = 1. / (6. * zeta_sketch_inter);
    let zeta;
    if sketch {
        zeta = zeta_sketch_inter;
    } else {
        zeta = zeta_inter;
    }
    let zeta = zeta * w;
    let zeta = zeta * 1000.;
    dbg!(zeta, c, k, n, s, alpha);

    let dump = true;
    if reference.contains(".sketch") {
        //        let reader = BufReader::new(File::open(&reference).expect(&reference));
        //        let res: Result<Sketch, _> = bincode::deserialize_from(reader);
        //        ref_sketch = res.unwrap();
    } else {
        std::mem::swap(
            &mut file_io::fastx_to_sketches(
                &vec![matches.value_of("reference").unwrap().to_string()],
                c,
                k,
            )[0],
            &mut ref_sketch,
        );
        if dump {
            let mut file_bin =
                BufWriter::new(File::create(format!("{}.sketch", &reference)).unwrap());
            bincode::serialize_into(&mut file_bin, &ref_sketch).unwrap();
        }
    }
    dbg!("done sketching reference");
    let read_sketch = file_io::fastx_to_multiple_sketch_rewrite(
        &vec![matches.value_of("query").unwrap().to_string()],
        c,
        k,
    );
    dbg!("done sketching query");
    let align_times: Mutex<Vec<_>> = Mutex::new(vec![]);

    (0..read_sketch.len())
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|i| {
            let query_sketch = &read_sketch[i];
            let times = align(&ref_sketch, query_sketch, k, w, zeta, debug);

            if !times.is_none() {
                let ct = times.unwrap().0;
                let et = times.unwrap().1;
                let covered = times.unwrap().2;
                align_times.lock().unwrap().push((
                    et,
                    ct,
                    query_sketch.contig_lengths[0],
                    covered,
                    query_sketch.contigs[0].clone(),
                ));
            }
        });

    let align_times = align_times.into_inner().unwrap();
    for align_time in align_times {
        println!(
            "{}\t{}\t{}\t{}\t{}",
            align_time.0, align_time.1, align_time.2, align_time.3, align_time.4
        );
    }
}

fn align(
    ref_sketch: &Sketch,
    query_sketch: &Sketch,
    k: usize,
    w: f64,
    zeta: f64,
    debug: bool,
) -> Option<(f64, f64, f64)> {
    let mut all_anchors = vec![];
    let now = Instant::now();
    let mut contig_count = vec![0; ref_sketch.contigs.len()];
    let kmer_seeds_query = query_sketch.kmer_seeds_k.as_ref().unwrap();
    let kmer_seeds_ref = ref_sketch.kmer_seeds_k.as_ref().unwrap();
    for (canon_kmer, query_pos) in kmer_seeds_query.iter() {
        let contains = kmer_seeds_ref.contains_key(canon_kmer);
        if contains {
            let ref_pos = &kmer_seeds_ref[canon_kmer];

            if ref_pos.len() > ref_sketch.repetitive_kmers {
                continue;
            }

            for qpos in query_pos {
                for rpos in ref_pos {
                    all_anchors.push(Anchor::new(
                        &(rpos.pos, rpos.contig_index),
                        &(qpos.pos, qpos.contig_index),
                        rpos.phase,
                        qpos.phase,
                        rpos.canonical != qpos.canonical,
                    ));
                }
            }
        }
    }

    if all_anchors.len() / query_sketch.contig_lengths[0] as usize > 10 {
        //Repetittive
        return Some((-3., -3., 0.));
    }
    let mut used_query_positions = FxHashMap::default();
    if all_anchors.len() > 0{
        for anchor in all_anchors.iter(){
            if used_query_positions.contains_key(&anchor.query_pos){
                used_query_positions.insert(anchor.query_pos,None);
            }
            else{
                used_query_positions.insert(anchor.query_pos, Some(anchor));
            }
        }
    }
    if all_anchors.len() > 0 {
        let mut rev_count = 0;
        for anchor_opt in used_query_positions.values(){
            if let Some(anchor) = anchor_opt{
                if anchor.reverse_match {
                    rev_count += 1;
                } else {
                    contig_count[anchor.ref_contig as usize] += 1;
                }
            }
        }
        if debug{
            dbg!(&contig_count);
        }
        //IMplementation isn't reverse-complement aware.
        if rev_count > all_anchors.len() / 2 {
            //            dbg!(rev_count,all_anchors.len(), kmer_seeds_query.len());
            return None;
        }
        let max_ind = seeding::position_max(&contig_count).unwrap();
        let ref_s = &ref_sketch.contig_seq[max_ind];
        let query_s = &query_sketch.contig_seq[0];
        let mut anchors = vec![];
        for anchor in all_anchors.iter() {
            if anchor.ref_contig as usize == max_ind && !anchor.reverse_match {
                anchors.push(anchor);
            }
        }
        if anchors.len() == 0 {
            return None;
        }
        anchors.sort();
        //        dbg!(&anchors);

        if debug {
            println!("Best contig is {}", &ref_sketch.contigs[max_ind]);
            println!("Number of anchors is {}", anchors.len());
            println!("Anchor finding time {}", now.elapsed().as_secs_f32());
        }

        let mut f = vec![w];
        let mut pointer_array = vec![0; anchors.len()];
        let mut avl_tree: SearchTree<[usize; 2]> = SearchTree::new();

        let now = Instant::now();
        for (i, anchor) in anchors.iter().enumerate() {
            avl_tree.insert([anchor.ref_pos as usize, i]);
        }
        avl_tree.update_query_info(
            [anchors[0].ref_pos as usize, 0],
            w + zeta * (anchors[0].ref_pos + anchors[0].query_pos) as f64,
            0,
            //don't think we need offset
            //            100 * anchors[0].query_pos as usize + offset,
            anchors[0].query_pos as usize,
            //            100 * anchors[0].ref_pos as usize + offset,
            anchors[0].ref_pos as usize,
        );

        for i in 1..anchors.len() {
            let best_f_i;
            let best_j;

            let (best_score, best_id) = avl_tree.mrq(
                [0, 0],
                [anchors[i].ref_pos as usize, i],
                anchors[i].query_pos as usize,
                anchors[i].ref_pos as usize,
            );
            if best_score == i64::MIN {
                best_f_i = w;
                best_j = i;
            } else {
                best_j = best_id;
                best_f_i = best_score as f64
                    - zeta * (anchors[i].query_pos + anchors[i].ref_pos) as f64
                    + w;
            }
            //                        if best_f_i < 0.0 {
            //                            best_f_i = 0.0;
            //                            best_j = i;
            //                        }
            f.push(best_f_i);
            avl_tree.update_query_info(
                [anchors[i].ref_pos as usize, i],
                best_f_i + zeta * (anchors[i].ref_pos + anchors[i].query_pos) as f64,
                i,
                anchors[i].query_pos as usize,
                anchors[i].ref_pos as usize,
            );

            if best_j != usize::MAX {
                pointer_array[i] = best_j;
            }
        }
        let chain_time = now.elapsed().as_secs_f64();
        let mut vec: Vec<_> = f.iter().enumerate().collect();
        vec.sort_by(|(_, v0), (_, v1)| v1.partial_cmp(v0).unwrap());
        let mut curr_i = vec[0].0;
        let mut prev_i = pointer_array[curr_i];
        let mut best_chain = vec![];
        while curr_i != prev_i {
            best_chain.push(anchors[curr_i]);
            curr_i = prev_i;
            prev_i = pointer_array[curr_i];
        }
        best_chain.push(anchors[curr_i]);
        let covered;
        if debug {
//            dbg!(&best_chain);
        }
        if best_chain.len() > 1 {
            let range = best_chain[0].ref_pos - best_chain[best_chain.len() - 1].ref_pos;
            let range_q = best_chain[0].query_pos - best_chain[best_chain.len() - 1].query_pos;
            covered = f64::min(1., range_q as f64 / query_sketch.contig_lengths[0] as f64);
            //            dbg!(range_q as f32 / query_sketch.contig_lengths[0] as f32, query_sketch.contig_lengths[0], &query_sketch.contigs[0]);
            if (range_q as f32) < 0.9 * query_sketch.contig_lengths[0] as f32 {
                return Some((-1., -1., covered));
            }
        } else {
            return None;
        }
        //            for anchors in best_chain.iter() {
        //                if anchors.0 != anchors.1 + start_ind {
        //                    if break_start == false {
        //                        break_start = true;
        //                        break_a = cmp::min(anchors.0, anchors.1);
        //                        break_b = cmp::max(anchors.0, anchors.1);
        //                    } else {
        //                        let min = cmp::min(anchors.0, anchors.1);
        //                        let max = cmp::max(anchors.0, anchors.1);
        //                        if min < break_a {
        //                            break_a = min;
        //                        }
        //                        if max > break_b {
        //                            break_b = max;
        //                        }
        //                    }
        //                } else {
        //                    if break_start {
        //                        break_start = false;
        //                        break_length += break_b - break_a;
        //                    }
        //                }
        //            }
        //            let mut rec = recov.lock().unwrap();
        //            if break_length > 0 {
        //                if print_all_debug {
        //                    dbg!(break_length);
        //                }
        //            }
        //            if range > break_length {
        //                *rec += (range - break_length) as f64;
        //            }
        //        }

        let mut gap_intervals = vec![];
        let mut curr_q_end = 0;
        let mut curr_r_end = 0;
        for (iter, anchor) in best_chain.iter().rev().enumerate() {
            if iter > 0 {
                let rgap;
                let qgap;
                if anchor.query_pos > curr_q_end + 1 && anchor.ref_pos > curr_r_end + 1 {
                    qgap = (curr_q_end, anchor.query_pos - 1);
                    rgap = (curr_r_end, anchor.ref_pos - 1);
                    gap_intervals.push((rgap, qgap));
                } else {
                }
            }
            curr_q_end = anchor.query_pos + k as u32 - 1;
            curr_r_end = anchor.ref_pos + k as u32 - 1;
        }

        //        let mut penalties = AffinePenalties {
        //            match_: 0,
        //            mismatch: 4,
        //            gap_opening: 6,
        //            gap_extension: 2,
        //        };

        //        dbg!(&gap_intervals);
        //        let alloc = MMAllocator::new(BUFFER_SIZE_8M as u64);
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let now = Instant::now();
        let mut aligner = Aligner::new(-5, -1, &score);
        let use_wfa = false;
        for gap in gap_intervals.iter() {
            let ref_slice = &ref_s[gap.0 .0 as usize..gap.0 .1 as usize];
            let query_slice = &query_s[gap.1 .0 as usize..gap.1 .1 as usize];
            //We get mistakes very rarely, but these create big outliers.
            if ref_slice.len() > 10000 || query_slice.len() > 10000 {
                //                dbg!(&gap_intervals);
                return Some((-2., -2., covered));
            }
            if !use_wfa {
                //rust bio only prints cigars for non-global alignment. uncomment
                //to verify that the alignment is actually correct.
                //                let alignment = aligner.semiglobal(&ref_slice, query_slice);
                //                dbg!(alignment.cigar(false));
                let alignment = aligner.global(&ref_slice, query_slice);
            } else {
                //                let mut wavefronts = AffineWavefronts::new_complete(
                //                    ref_slice.len(),
                //                    query_slice.len(),
                //                    &mut penalties,
                //                    &alloc,
                //                );
                //                wavefronts.align(ref_slice, query_slice).unwrap();
            }
        }

        let extend_time = now.elapsed().as_secs_f64();
        return Some((chain_time, extend_time, covered));
    }
    return None;
}
