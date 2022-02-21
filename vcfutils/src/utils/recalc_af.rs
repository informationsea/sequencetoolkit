// Recalculate AF/AC/AN

use crate::error::{VCFUtilsError, VCFUtilsErrorKind};
use nom::{
    self, branch::alt, bytes::complete::tag, bytes::complete::take_while1, character::is_digit,
    eof, multi::separated_list, named, sequence::tuple,
};
use std::str;
use vcf::{U8Vec, VCFRecord};

named!(eof_parser, eof!());

pub struct AlleleCount {
    pub allele_count: Vec<u64>,
    pub allele_number: u64,
    pub genotype_count: Vec<u64>,
    pub alt_hom_count: Vec<u64>,
    pub ploidy: usize,
    cache: Vec<usize>,
}

fn is_digit_or_dot(chr: u8) -> bool {
    is_digit(chr) || chr == b'.'
}

impl AlleleCount {
    pub fn new(alt_count: usize, ploidy: usize) -> AlleleCount {
        let mut allele_count = Vec::new();
        let mut alt_hom_count = Vec::new();
        for _ in 0..(alt_count + 1) {
            allele_count.push(0);
            alt_hom_count.push(0);
        }
        let genotype_count = vec![0; number_of_genotype(alt_count, ploidy)];

        AlleleCount {
            allele_count,
            allele_number: 0,
            genotype_count,
            alt_hom_count,
            ploidy,
            cache: Vec::new(),
        }
    }

    pub fn add_genotype(&mut self, genotype: &[u8]) -> Result<(), VCFUtilsError> {
        match tuple((
            separated_list(alt((tag("|"), tag("/"))), take_while1(is_digit_or_dot)),
            eof_parser,
        ))(genotype)
        {
            Ok((_, (parsed_genotype, _))) => {
                self.cache.clear();
                for one_allele in &parsed_genotype {
                    if one_allele != b"." {
                        let allele: usize = str::from_utf8(one_allele)?.parse::<usize>()?;
                        self.cache.push(allele);
                        if let Some(x) = self.allele_count.get_mut(allele) {
                            *x += 1;
                            self.allele_number += 1;
                        }
                    }
                }
                self.cache.sort_unstable();
                if self.cache.len() == self.ploidy {
                    if let Some(x) = self.genotype_count.get_mut(genotype_index(&self.cache)) {
                        *x += 1;
                    }
                    if let Some(first_allele) = self.cache.get(0) {
                        if self.cache.iter().all(|x| x == first_allele) {
                            if let Some(x) = self.alt_hom_count.get_mut(*first_allele) {
                                *x += 1;
                            }
                        }
                    }
                }

                Ok(())
            }
            Err(_) => Err(VCFUtilsErrorKind::GenotypeParseError.into()),
        }
    }

    pub fn add_record<'a, I: Iterator<Item = &'a U8Vec>>(
        &mut self,
        record: &VCFRecord,
        ploidy: usize,
        samples: I,
    ) -> Result<(), VCFUtilsError> {
        self.clear(record.alternative.len(), ploidy);
        for one_sample in samples {
            record
                .genotype(one_sample, b"GT")
                .map(|x| x.get(0))
                .flatten()
                .map(|x| self.add_genotype(x));
        }
        Ok(())
    }

    pub fn clear(&mut self, new_alt_allele_count: usize, new_ploidy: usize) {
        self.allele_count.clear();
        self.alt_hom_count.clear();
        for _ in 0..(new_alt_allele_count + 1) {
            self.allele_count.push(0);
            self.alt_hom_count.push(0);
        }
        self.genotype_count.clear();
        self.genotype_count
            .resize(number_of_genotype(new_alt_allele_count, new_ploidy), 0);
        self.allele_number = 0;
    }
}

pub fn number_of_genotype(alt_allele_num: usize, ploidy: usize) -> usize {
    1 + (0..ploidy)
        .map(|m| combination(m + 1, alt_allele_num + m))
        .sum::<usize>()
}

/// alleles should be sorted
pub fn genotype_index(alleles: &[usize]) -> usize {
    (0..alleles.len())
        .map(|m| combination(m + 1, alleles[m] + m))
        .sum()
}

pub fn combination(k: usize, n: usize) -> usize {
    if n < k {
        0
    } else if k == 0 {
        1
    } else if (n - k) < k {
        combination(n - k, n)
    } else {
        combination(k - 1, n - 1) * n / k
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::error::VCFUtilsError;
    use std::io::BufReader;
    use vcf::VCFReader;

    #[test]
    fn test_combination() {
        assert_eq!(combination(1, 1), 1);
        assert_eq!(combination(0, 1), 1);
        assert_eq!(combination(1, 1), 1);
        assert_eq!(combination(1, 3), 3);
        assert_eq!(combination(2, 3), 3);
        assert_eq!(combination(1, 4), 4);
        assert_eq!(combination(2, 4), 6);
    }

    #[test]
    fn test_number_of_genotype() {
        assert_eq!(number_of_genotype(1, 2), 3);
        assert_eq!(number_of_genotype(2, 2), 6);
        assert_eq!(number_of_genotype(3, 2), 10);
        assert_eq!(number_of_genotype(2, 3), 10);
    }

    #[test]
    fn test_genotype_index() {
        assert_eq!(genotype_index(&[0, 0]), 0);
        assert_eq!(genotype_index(&[0, 1]), 1);
        assert_eq!(genotype_index(&[1, 1]), 2);
        assert_eq!(genotype_index(&[0, 2]), 3);
        assert_eq!(genotype_index(&[1, 2]), 4);
        assert_eq!(genotype_index(&[2, 2]), 5);

        assert_eq!(genotype_index(&[0, 0, 0]), 0);
        assert_eq!(genotype_index(&[0, 0, 1]), 1);
        assert_eq!(genotype_index(&[0, 1, 1]), 2);
        assert_eq!(genotype_index(&[1, 1, 1]), 3);
        assert_eq!(genotype_index(&[0, 0, 2]), 4);
        assert_eq!(genotype_index(&[0, 1, 2]), 5);
        assert_eq!(genotype_index(&[1, 1, 2]), 6);
        assert_eq!(genotype_index(&[0, 2, 2]), 7);
        assert_eq!(genotype_index(&[1, 2, 2]), 8);
        assert_eq!(genotype_index(&[2, 2, 2]), 9);
    }

    #[test]
    fn test_from_record() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../testfiles/simple1.vcf");
        let mut data_reader = BufReader::new(&vcf_data[..]);
        let mut vcf_reader = VCFReader::new(&mut data_reader)?;
        let mut record = VCFRecord::new(vcf_reader.header().clone());
        let mut allele_count = AlleleCount::new(1, 2);

        vcf_reader.next_record(&mut record)?;
        allele_count.add_record(&record, 2, vcf_reader.header().samples().iter())?;
        assert_eq!(allele_count.allele_count, vec![0, 4]);
        assert_eq!(allele_count.allele_number, 4);
        assert_eq!(allele_count.genotype_count, vec![0, 0, 2]);
        assert_eq!(allele_count.alt_hom_count, vec![0, 2]);

        vcf_reader.next_record(&mut record)?;
        allele_count.add_record(&record, 2, vcf_reader.header().samples().iter())?;
        assert_eq!(allele_count.allele_count, vec![3, 1]);
        assert_eq!(allele_count.allele_number, 4);
        assert_eq!(allele_count.genotype_count, vec![1, 1, 0]);
        assert_eq!(allele_count.alt_hom_count, vec![1, 0]);

        vcf_reader.next_record(&mut record)?;
        allele_count.add_record(&record, 2, vcf_reader.header().samples().iter())?;
        assert_eq!(allele_count.allele_count, vec![0, 2, 2]);
        assert_eq!(allele_count.allele_number, 4);
        assert_eq!(allele_count.genotype_count, vec![0, 0, 0, 0, 2, 0]);
        assert_eq!(allele_count.alt_hom_count, vec![0, 0, 0]);

        Ok(())
    }

    #[test]
    fn test_allele_count() -> Result<(), VCFUtilsError> {
        let mut allele_count = AlleleCount::new(2, 2);
        allele_count.add_genotype(b"0/0")?;
        allele_count.add_genotype(b"0/0")?;
        allele_count.add_genotype(b"0/1")?;
        allele_count.add_genotype(b"1/2")?;
        assert_eq!(allele_count.allele_count, vec![5, 2, 1]);
        assert_eq!(allele_count.genotype_count, vec![2, 1, 0, 0, 1, 0]);
        assert_eq!(allele_count.alt_hom_count, vec![2, 0, 0]);

        Ok(())
    }

    #[test]
    fn test_allele_count2() -> Result<(), VCFUtilsError> {
        let mut allele_count = AlleleCount::new(2, 2);
        allele_count.add_genotype(b"0/0")?;
        allele_count.add_genotype(b"2/2")?;
        allele_count.add_genotype(b"0/1")?;
        allele_count.add_genotype(b"1/1")?;
        assert_eq!(allele_count.allele_count, vec![3, 3, 2]);
        assert_eq!(allele_count.genotype_count, vec![1, 1, 1, 0, 0, 1]);
        assert_eq!(allele_count.alt_hom_count, vec![1, 1, 1]);

        Ok(())
    }

    #[test]
    fn test_allele_count3() -> Result<(), VCFUtilsError> {
        let mut allele_count = AlleleCount::new(2, 2);
        allele_count.add_genotype(b"0/0")?;
        allele_count.add_genotype(b"./.")?;
        allele_count.add_genotype(b"0/1")?;
        allele_count.add_genotype(b"1/1")?;
        assert_eq!(allele_count.allele_count, vec![3, 3, 0]);
        assert_eq!(allele_count.genotype_count, vec![1, 1, 1, 0, 0, 0]);
        assert_eq!(allele_count.alt_hom_count, vec![1, 1, 0]);

        Ok(())
    }
}
