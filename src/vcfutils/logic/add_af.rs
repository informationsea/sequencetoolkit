use crate::vcfutils::error::VCFUtilsError;
use crate::vcfutils::utils::recalc_af::AlleleCount;
use itertools::concat;
use std::collections::{HashMap, HashSet};
use std::io::prelude::*;
use std::str;
use vcf::{self, U8Vec};

pub fn add_af<R: BufRead, W: Write>(
    reader: &mut vcf::VCFReader<R>,
    writer: W,
    category_to_sample: &CategoryToSamples,
) -> Result<(), VCFUtilsError> {
    let mut info_keys_to_samples: Vec<(AlleleCountInfoKeys, HashSet<U8Vec>)> = category_to_sample
        .iter()
        .map(|(k, v)| {
            (
                AlleleCountInfoKeys {
                    ac: concat(vec![b"AC_".to_vec(), k.to_vec()]),
                    an: concat(vec![b"AN_".to_vec(), k.to_vec()]),
                    af: concat(vec![b"AF_".to_vec(), k.to_vec()]),
                    genotype_count: concat(vec![b"GenotypeCount_".to_vec(), k.to_vec()]),
                    alt_hom_count: concat(vec![b"nhomalt_".to_vec(), k.to_vec()]),
                },
                v.clone(),
            )
        })
        .collect();
    info_keys_to_samples.insert(
        0,
        (
            AlleleCountInfoKeys {
                ac: b"AC".to_vec(),
                an: b"AN".to_vec(),
                af: b"AF".to_vec(),
                genotype_count: b"GenotypeCount".to_vec(),
                alt_hom_count: b"nhomalt".to_vec(),
            },
            reader.header().samples().iter().cloned().collect(),
        ),
    );

    let mut header_items = reader.header().items().to_vec();
    for (one_info_keys, _) in info_keys_to_samples.iter() {
        if reader.header().info(&one_info_keys.ac).is_none() {
            header_items.push(vcf::VCFHeaderLine::from_bytes(
                &format!(
                    "##INFO=<ID={},Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n",
                    str::from_utf8(&one_info_keys.ac).unwrap()
                ).into_bytes(),
            0).unwrap());
        }
        if reader.header().info(&one_info_keys.an).is_none() {
            header_items.push(vcf::VCFHeaderLine::from_bytes(
                &format!(
                    "##INFO=<ID={},Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n",
                    str::from_utf8(&one_info_keys.an).unwrap()
                ).into_bytes(),
            0).unwrap());
        }
        if reader.header().info(&one_info_keys.af).is_none() {
            header_items.push(vcf::VCFHeaderLine::from_bytes(
                &format!(
                    "##INFO=<ID={},Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n",
                    str::from_utf8(&one_info_keys.af).unwrap()
                ).into_bytes(),
            0).unwrap());
        }
        if reader
            .header()
            .info(&one_info_keys.genotype_count)
            .is_none()
        {
            header_items.push(vcf::VCFHeaderLine::from_bytes(
                &format!(
                    "##INFO=<ID={},Number=G,Type=Integer,Description=\"Genotype count in called genotypes. Refer to Genotype Ordering section in VCF v4.3 specification for ordering of genotypes.\">\n",
                    str::from_utf8(&one_info_keys.genotype_count).unwrap()
                ).into_bytes(),
            0).unwrap());
        }
        if reader.header().info(&one_info_keys.alt_hom_count).is_none() {
            header_items.push(vcf::VCFHeaderLine::from_bytes(
                &format!(
                    "##INFO=<ID={},Number=A,Type=Integer,Description=\"Count of homozygous individuals\">\n",
                    str::from_utf8(&one_info_keys.alt_hom_count).unwrap()
                ).into_bytes(),
            0).unwrap());
        }
    }
    let new_header = vcf::VCFHeader::new(header_items, reader.header().samples().to_vec());

    let mut record = vcf::VCFRecord::new(reader.header().clone());
    let mut allele_count = AlleleCount::new(1, 2);
    let mut vcf_writer = vcf::VCFWriter::new(writer, &new_header)?;

    while reader.next_record(&mut record)? {
        for (k, v) in info_keys_to_samples.iter() {
            add_af_to_record(&mut record, k, v, &mut allele_count)?;
        }
        vcf_writer.write_record(&record)?;
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct AlleleCountInfoKeys {
    pub ac: U8Vec,
    pub an: U8Vec,
    pub af: U8Vec,
    pub genotype_count: U8Vec,
    pub alt_hom_count: U8Vec,
}

pub fn add_af_to_record<S>(
    record: &mut vcf::VCFRecord,
    info_keys: &AlleleCountInfoKeys,
    samples: &HashSet<U8Vec, S>,
    cache: &mut AlleleCount,
) -> Result<(), VCFUtilsError> {
    cache.add_record(record, 2, samples.iter())?;
    record.insert_info(
        &info_keys.an,
        vec![format!("{}", cache.allele_number).into_bytes()],
    );
    record.insert_info(
        &info_keys.ac,
        cache
            .allele_count
            .iter()
            .skip(1)
            .map(|x| format!("{}", x).into_bytes())
            .collect(),
    );
    record.insert_info(
        &info_keys.af,
        cache
            .allele_count
            .iter()
            .skip(1)
            .map(|x| format!("{:.4}", (*x as f64) / (cache.allele_number as f64)).into_bytes())
            .collect(),
    );
    record.insert_info(
        &info_keys.genotype_count,
        cache
            .genotype_count
            .iter()
            .map(|x| format!("{}", x).into_bytes())
            .collect(),
    );
    record.insert_info(
        &info_keys.alt_hom_count,
        cache
            .alt_hom_count
            .iter()
            .skip(1)
            .map(|x| format!("{}", x).into_bytes())
            .collect(),
    );
    Ok(())
}

type CategoryToSamples = HashMap<U8Vec, HashSet<U8Vec>>;

pub fn load_category_mapping<R: Read, S: std::hash::BuildHasher>(
    reader: &mut csv::Reader<R>,
    id_keys: HashSet<U8Vec, S>,
    category_keys: HashSet<U8Vec, S>,
) -> Result<CategoryToSamples, VCFUtilsError> {
    let bytes_header = reader.byte_headers()?;
    let header: Vec<U8Vec> = bytes_header.iter().map(|x| x.to_vec()).collect();
    let id_index: HashSet<usize> = header
        .iter()
        .enumerate()
        .filter(|x| id_keys.contains(x.1))
        .map(|x| x.0)
        .collect();
    let category_index: HashSet<usize> = header
        .iter()
        .enumerate()
        .filter(|x| category_keys.contains(x.1))
        .map(|x| x.0)
        .collect();
    let mut category_mapping: CategoryToSamples = HashMap::new();

    let mut record = csv::ByteRecord::new();
    while reader.read_byte_record(&mut record)? {
        for one_id in id_index.iter().filter_map(|x| record.get(*x)) {
            for one_category in category_index.iter().filter_map(|x| record.get(*x)) {
                if let Some(x) = category_mapping.get_mut(one_category) {
                    x.insert(one_id.to_vec());
                } else {
                    let mut s = HashSet::new();
                    s.insert(one_id.to_vec());
                    category_mapping.insert(one_category.to_vec(), s);
                }
            }
        }
    }

    Ok(category_mapping)
}

#[cfg(test)]
mod test {
    use super::*;
    use std::collections::hash_map::RandomState;

    #[allow(clippy::cognitive_complexity)]
    #[test]
    fn test_add_af() -> Result<(), VCFUtilsError> {
        let vcf_data = include_bytes!("../../../testfiles/vcfutils/1kGP-subset.vcf");
        let mut vcf_reader = vcf::VCFReader::new(&vcf_data[..])?;
        let category_mapping = load_category_mapping::<_, RandomState>(
            &mut csv::Reader::from_reader(
                &include_bytes!("../../../testfiles/vcfutils/1kGP-category.csv")[..],
            ),
            [b"id".to_vec()].iter().cloned().collect(),
            [b"category1".to_vec(), b"category2".to_vec()]
                .iter()
                .cloned()
                .collect(),
        )?;
        let mut data = Vec::<u8>::new();
        add_af(&mut vcf_reader, &mut data, &category_mapping)?;
        std::fs::File::create("target/add-af.vcf")?.write_all(&data)?;

        let mut vcf_reader2 = vcf::VCFReader::new(&data[..])?;
        let info_set = vcf_reader2
            .header()
            .info_list()
            .cloned()
            .collect::<HashSet<_>>();

        // check info
        assert_eq!(true, info_set.contains(&b"AC".to_vec()));
        assert_eq!(true, info_set.contains(&b"AN".to_vec()));
        assert_eq!(true, info_set.contains(&b"AF".to_vec()));
        assert_eq!(true, info_set.contains(&b"GenotypeCount".to_vec()));
        assert_eq!(true, info_set.contains(&b"AC_P1".to_vec()));
        assert_eq!(true, info_set.contains(&b"AN_P1".to_vec()));
        assert_eq!(true, info_set.contains(&b"AF_P1".to_vec()));
        assert_eq!(true, info_set.contains(&b"GenotypeCount_P1".to_vec()));
        assert_eq!(true, info_set.contains(&b"AC_P2".to_vec()));
        assert_eq!(true, info_set.contains(&b"AN_P2".to_vec()));
        assert_eq!(true, info_set.contains(&b"AF_P2".to_vec()));
        assert_eq!(true, info_set.contains(&b"GenotypeCount_P2".to_vec()));
        assert_eq!(true, info_set.contains(&b"AC_X1".to_vec()));
        assert_eq!(true, info_set.contains(&b"AN_X1".to_vec()));
        assert_eq!(true, info_set.contains(&b"AF_X1".to_vec()));
        assert_eq!(true, info_set.contains(&b"GenotypeCount_X1".to_vec()));
        assert_eq!(true, info_set.contains(&b"AC_X2".to_vec()));
        assert_eq!(true, info_set.contains(&b"AN_X2".to_vec()));
        assert_eq!(true, info_set.contains(&b"AF_X2".to_vec()));
        assert_eq!(true, info_set.contains(&b"GenotypeCount_X2".to_vec()));

        let mut vcf_record = vcf::VCFRecord::new(vcf_reader2.header().clone());
        assert_eq!(true, vcf_reader2.next_record(&mut vcf_record)?);
        assert_eq!(vcf_record.info(b"AC_P1"), Some(&vec![b"1".to_vec()]));
        assert_eq!(vcf_record.info(b"AC_P2"), Some(&vec![b"0".to_vec()]));
        assert_eq!(
            vcf_record.info(b"GenotypeCount_P1"),
            Some(&vec![b"1".to_vec(), b"1".to_vec(), b"0".to_vec()])
        );
        assert_eq!(
            vcf_record.info(b"GenotypeCount_P2"),
            Some(&vec![b"3".to_vec(), b"0".to_vec(), b"0".to_vec()])
        );

        assert_eq!(true, vcf_reader2.next_record(&mut vcf_record)?);
        assert_eq!(vcf_record.info(b"AC_P1"), Some(&vec![b"2".to_vec()]));
        assert_eq!(vcf_record.info(b"AC_P2"), Some(&vec![b"3".to_vec()]));
        assert_eq!(
            vcf_record.info(b"GenotypeCount_P1"),
            Some(&vec![b"0".to_vec(), b"2".to_vec(), b"0".to_vec()])
        );
        assert_eq!(
            vcf_record.info(b"GenotypeCount_P2"),
            Some(&vec![b"1".to_vec(), b"1".to_vec(), b"1".to_vec()])
        );

        assert_eq!(true, vcf_reader2.next_record(&mut vcf_record)?);
        assert_eq!(vcf_record.info(b"AC_P1"), Some(&vec![b"2".to_vec()]));
        assert_eq!(vcf_record.info(b"AC_P2"), Some(&vec![b"1".to_vec()]));
        assert_eq!(
            vcf_record.info(b"GenotypeCount_P1"),
            Some(&vec![b"0".to_vec(), b"2".to_vec(), b"0".to_vec()])
        );
        assert_eq!(
            vcf_record.info(b"GenotypeCount_P2"),
            Some(&vec![b"2".to_vec(), b"1".to_vec(), b"0".to_vec()])
        );

        Ok(())
    }

    #[test]
    fn test_load_category_mapping() {
        let data = include_bytes!("../../../testfiles/vcfutils/category.csv");
        let mut csv_reader = csv::Reader::from_reader(&data[..]);
        let result = load_category_mapping::<_, RandomState>(
            &mut csv_reader,
            [b"id1".to_vec(), b"id2".to_vec()].iter().cloned().collect(),
            [b"sex".to_vec(), b"protocol".to_vec()]
                .iter()
                .cloned()
                .collect(),
        )
        .unwrap();
        assert_eq!(
            result,
            [
                (
                    b"P1".to_vec(),
                    [b"A".to_vec(), b"1".to_vec(), b"C".to_vec(), b"3".to_vec()]
                        .iter()
                        .cloned()
                        .collect()
                ),
                (
                    b"P2".to_vec(),
                    [b"B".to_vec(), b"2".to_vec(), b"D".to_vec(), b"4".to_vec()]
                        .iter()
                        .cloned()
                        .collect()
                ),
                (
                    b"F".to_vec(),
                    [b"A".to_vec(), b"1".to_vec(), b"B".to_vec(), b"2".to_vec()]
                        .iter()
                        .cloned()
                        .collect()
                ),
                (
                    b"M".to_vec(),
                    [b"C".to_vec(), b"3".to_vec(), b"D".to_vec(), b"4".to_vec()]
                        .iter()
                        .cloned()
                        .collect()
                )
            ]
            .iter()
            .cloned()
            .collect(),
        )
    }
}
