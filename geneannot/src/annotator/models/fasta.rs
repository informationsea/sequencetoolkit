use super::{Chromosome, Genome};
use bio::io::fasta::Index;

pub fn load_fasta(name: &str, index: &Index) -> Genome {
    let chromosomes: Vec<_> = index
        .sequences()
        .iter()
        .map(|x| Chromosome {
            name: x.name.to_string(),
            length: x.len,
        })
        .collect();

    Genome::new(name, &chromosomes[..])
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_load_fasta() {
        let index = Index::from_file(
            &"testfiles/genome/Influenza/GCF_000865085.1_ViralMultiSegProj15622_genomic.fna.fai",
        )
        .unwrap();
        assert_eq!(
            load_fasta("influenza", &index),
            Genome {
                name: "influenza".to_string(),
                chromosomes: vec![
                    Chromosome {
                        name: "NC_007373.1".to_string(),
                        length: 2341
                    },
                    Chromosome {
                        name: "NC_007372.1".to_string(),
                        length: 2341
                    },
                    Chromosome {
                        name: "NC_007371.1".to_string(),
                        length: 2233
                    },
                    Chromosome {
                        name: "NC_007366.1".to_string(),
                        length: 1762
                    },
                    Chromosome {
                        name: "NC_007369.1".to_string(),
                        length: 1566
                    },
                    Chromosome {
                        name: "NC_007368.1".to_string(),
                        length: 1467
                    },
                    Chromosome {
                        name: "NC_007367.1".to_string(),
                        length: 1027
                    },
                    Chromosome {
                        name: "NC_007370.1".to_string(),
                        length: 890
                    }
                ],
                name_to_index: [
                    ("NC_007373.1", 0),
                    ("NC_007372.1", 1),
                    ("NC_007371.1", 2),
                    ("NC_007366.1", 3),
                    ("NC_007369.1", 4),
                    ("NC_007368.1", 5),
                    ("NC_007367.1", 6),
                    ("NC_007370.1", 7)
                ]
                .iter()
                .map(|(k, v)| ((*k).to_string(), *v))
                .collect()
            }
        )
    }
}
